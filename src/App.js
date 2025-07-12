import React, { useState, useCallback } from 'react';
// import * as ReactDOM from 'react-dom'; // Removed as it's not directly used in App component

// Utility function for lstrip (moved from String.prototype)
const lstripString = (str, chars) => {
    if (!chars) {
        return str.trimStart(); // Use trimStart for modern JS equivalent of trimLeft
    }
    let pattern = new RegExp('^[' + chars + ']+');
    return str.replace(pattern, '');
};

function App() {
    const [fastaFile, setFastaFile] = useState(null);
    const [statusMessage, setStatusMessage] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const [resultFasta, setResultFasta] = useState('');
    const [error, setError] = useState('');
    const [variantCriteria, setVariantCriteria] = useState('pathogenicOnly'); // New state for variant selection
    const [showUtilityExample, setShowUtilityExample] = useState(false); // New state for toggling utility example

    // Function to extract UniProt IDs and sequences from a FASTA file
    const getSequencesAndIdsFromFasta = useCallback((fastaContent) => {
        const sequences = {};
        const uniprotIds = [];
        const canonicalHeaders = {};
        let currentId = null;
        let currentSeqLines = [];

        const lines = fastaContent.split('\n');
        for (const line of lines) {
            const trimmedLine = line.trim();
            if (trimmedLine.startsWith('>')) {
                if (currentId && currentSeqLines.length > 0) {
                    sequences[currentId] = currentSeqLines.join('');
                }

                // MODIFIED REGEX HERE to correctly extract base UniProt ID
                // This regex captures the main UniProt ID (e.g., P08246) and optionally ignores -VAR_XXXXXX suffixes.
                const match = trimmedLine.match(/^>[a-z]{2}\|([A-Z0-9_.-]+?)(?:-VAR_[0-9]+)?\|/);
                if (match) {
                    currentId = match[1]; // This will now correctly be P08246
                } else {
                    // Fallback for less standard headers, still try to get the base ID
                    const parts = trimmedLine.substring(1).split(' ');
                    if (parts.length > 0) {
                        let potentialId = parts[0];
                        // Try to extract the base ID from potentialId if it contains VAR_
                        const varMatch = potentialId.match(/^([A-Z0-9_.-]+?)(?:-VAR_[0-9]+)?$/);
                        if (varMatch) {
                            currentId = varMatch[1];
                        } else {
                            currentId = potentialId;
                        }
                    } else {
                        currentId = `UNKNOWN_ID_${Math.random().toString(36).substring(7)}`; // Generate random ID
                    }
                }

                if (currentId && !uniprotIds.includes(currentId)) {
                    uniprotIds.push(currentId);
                }
                // Store the complete original header line (including VAR_ if present in input)
                canonicalHeaders[currentId] = trimmedLine; 
                currentSeqLines = [];
            } else {
                currentSeqLines.push(trimmedLine);
            }
        }
        if (currentId && currentSeqLines.length > 0) {
            sequences[currentId] = currentSeqLines.join('');
        }
        return { uniprotIds, sequences, canonicalHeaders };
    }, []);

    // Function to process a protein and its variants
    const processProteinVariants = useCallback((uniprotId, canonicalSequence, canonicalHeaderLine, variantsData, currentVariantCriteria) => {
        const pathogenicVariantSequences = [];
        if (!variantsData || !variantsData.features) {
            return [];
        }

        for (const feature of variantsData.features) {
            // Condition to include the variant based on selected criteria
            let includeFeature = false;
            let isPathogenic = false; // Determine pathogenicity regardless of inclusion criteria for description/ftId

            if (feature.clinicalSignificances) {
                for (const cs of feature.clinicalSignificances) {
                    if (cs.type === "Pathogenic" || cs.type === "Likely pathogenic") {
                        isPathogenic = true;
                        break;
                    }
                }
            }

            if (currentVariantCriteria === 'pathogenicOnly') {
                if ((feature.type === "VARIANT" || feature.type === "Natural variant") && isPathogenic) {
                    includeFeature = true;
                }
            } else if (currentVariantCriteria === 'allWithDescription') {
                // Include if it's a VARIANT or Natural variant AND has a description
                // Now also check for clinicalSignificances or diseaseAssociations or association
                if ((feature.type === "VARIANT" || feature.type === "Natural variant") && 
                    (feature.description || (feature.clinicalSignificances && feature.clinicalSignificances.length > 0) || (feature.diseaseAssociations && feature.diseaseAssociations.length > 0) || (feature.association && feature.association.length > 0))) {
                    includeFeature = true;
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Included variant for 'All Variants with Descriptions': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Description: "${feature.description}"`);
                    // --- END DEBUG LOG ---
                } else {
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Skipped variant for 'All Variants with Descriptions': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Has Description: ${!!feature.description}, Has Clinical Significances: ${!!(feature.clinicalSignificances && feature.clinicalSignificances.length > 0)}, Has Disease Associations: ${!!(feature.diseaseAssociations && feature.diseaseAssociations.length > 0)}, Has Association: ${!!(feature.association && feature.association.length > 0)}`);
                    // --- END DEBUG LOG ---
                }
            } else if (currentVariantCriteria === 'allVariants') {
                // Include all 'VARIANT' or 'Natural variant' features
                if (feature.type === "VARIANT" || feature.type === "Natural variant") {
                    includeFeature = true;
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Included variant for 'All Variants': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Has Description: ${!!feature.description}`);
                    // --- END DEBUG ---
                }
            }

            if (includeFeature) {
                const beginPos = parseInt(feature.begin);
                const endPos = parseInt(feature.end);
                
                // Be more robust with wildType/mutatedType, falling back to alternativeSequence if direct fields are missing
                let wildType = feature.wildType || (feature.alternativeSequence && feature.alternativeSequence.originalSequence);
                let mutatedType = feature.mutatedType || (feature.alternativeSequence && feature.alternativeSequence.alternativeSequences && feature.alternativeSequence.alternativeSequences[0]);

                // --- NEW LOGIC: Handle '?' in mutatedType by replacing with 'X' ---
                if (mutatedType === '?') {
                    console.warn(`[WARNING] Replaced '?' mutatedType with 'X' for variant at position ${beginPos} of ${uniprotId}.`);
                    mutatedType = 'X';
                }
                // --- END NEW LOGIC ---

                let mutationDescription = "";
                if (wildType && mutatedType) {
                    mutationDescription = `p.${wildType}${beginPos}${mutatedType}`;
                } else if (wildType && !mutatedType && beginPos === endPos) {
                    mutationDescription = `p.del${wildType}${beginPos}`;
                } else if (wildType && !mutatedType && beginPos !== endPos) {
                    mutationDescription = `p.del${wildType}${beginPos}-${endPos}`;
                } else if (!wildType && mutatedType) { // Insertion
                    mutationDescription = `p.ins${mutatedType}${beginPos}`;
                }
                // Default mutation description if no specific type matches
                if (!mutationDescription) {
                    mutationDescription = `Variant at pos ${beginPos}`;
                }

                let modifiedSequence = ''; 
                let variantSequenceList = Array.from(canonicalSequence);
                let mutationAppliedSuccessfully = false;

                // --- MODIFIED MUTATION APPLICATION LOGIC ---
                if (beginPos - 1 < 0 || beginPos - 1 >= variantSequenceList.length) {
                    console.warn(`Warning: Mutation position ${beginPos} for ${uniprotId} is out of bounds for sequence length ${canonicalSequence.length}. Skipping mutation.`);
                } else {
                    const currentAminoAcidInFasta = variantSequenceList[beginPos - 1]; // Corrected variable name

                    if (wildType && mutatedType) { // Substitution
                        if (currentAminoAcidInFasta !== wildType) { // Corrected variable name
                            console.warn(`Warning: For ${uniprotId} at pos ${beginPos}, canonical AA in FASTA ('${currentAminoAcidInFasta}') does not match UniProt wildType ('${wildType}'). Applying mutation anyway.`);
                        }
                        variantSequenceList[beginPos - 1] = mutatedType;
                        mutationAppliedSuccessfully = true;
                    } else if (wildType && !mutatedType) { // Deletion
                        if (endPos <= variantSequenceList.length && beginPos <= endPos) {
                            const expectedDeletedInFasta = canonicalSequence.substring(beginPos - 1, endPos);
                            if (expectedDeletedInFasta !== wildType) {
                                console.warn(`Warning: For ${uniprotId} at pos ${beginPos}-${endPos}, canonical sequence in FASTA ('${expectedDeletedInFasta}') does not match UniProt wildType ('${wildType}'). Applying deletion anyway.`);
                            }
                            variantSequenceList.splice(beginPos - 1, endPos - (beginPos - 1));
                            mutationAppliedSuccessfully = true;
                        } else {
                            console.warn(`Warning: Deletion position ${beginPos}-${endPos} for ${uniprotId} is out of bounds or invalid. Skipping mutation.`);
                        }
                    } else if (!wildType && mutatedType) { // Insertion
                        if (beginPos - 1 <= variantSequenceList.length) {
                            variantSequenceList.splice(beginPos - 1, 0, mutatedType);
                            mutationAppliedSuccessfully = true;
                        } else {
                            console.warn(`Warning: Insertion position ${beginPos} for ${uniprotId} is out of bounds. Skipping mutation.`);
                        }
                    }
                }
                // --- END MODIFIED MUTATION APPLICATION LOGIC ---

                if (mutationAppliedSuccessfully) {
                    modifiedSequence = variantSequenceList.join('');
                    
                    const genomicLocInfo = feature.genomicLocation || [''];
                    const firstGenomicLoc = genomicLocInfo[0] || '';

                    // Use lstripString to remove the leading '>'
                    const cleanedOriginalHeader = lstripString(canonicalHeaderLine.trim(), '>');
                    const headerParts = cleanedOriginalHeader.split('|', 3); // Split into at most 3 parts: sp, ID, rest

                    let finalAccessionPart = uniprotId; // Default to base UniProt ID

                    if (headerParts.length >= 2) {
                        const originalAccessionInHeader = headerParts[1]; // e.g., P08246 or P08246-VAR_000100

                        // Check if the original accession already contains "-VAR_"
                        if (!originalAccessionInHeader.includes('-VAR_')) {
                            // If it doesn't have -VAR_, append the mutation description (without "p.")
                            finalAccessionPart = `${originalAccessionInHeader}-${mutationDescription.replace('p.', '')}`;
                        } else {
                            // If it already has -VAR_, keep it as is
                            finalAccessionPart = originalAccessionInHeader;
                        }

                        // Reconstruct the modified header part with the new ID
                        headerParts[1] = finalAccessionPart;
                        modifiedCanonicalHeaderPart = headerParts.join('|');
                    } else {
                        // Fallback for very simple header formats (e.g., just ">ID")
                        // In this case, always add the mutation to ensure uniqueness
                        modifiedCanonicalHeaderPart = `${uniprotId}-${mutationDescription.replace('p.', '')}`;
                        if (headerParts.length > 0) {
                            modifiedCanonicalHeaderPart = `${headerParts[0]}|${modifiedCanonicalHeaderPart}`;
                        }
                    }

                    // Determine the appropriate header tag based on pathogenicity
                    const headerTag = isPathogenic ? 'PATHOGENIC_VARIANT' : 'VARIANT_INFO'; 

                    // --- MODIFIED LOGIC FOR FEATURE DESCRIPTION PART ---
                    let featureDescriptionPart = '';
                    console.log(`[DEBUG] Building description for ${uniprotId} variant at pos ${beginPos}.`);
                    console.log(`[DEBUG]   feature.diseaseAssociations:`, feature.diseaseAssociations);
                    console.log(`[DEBUG]   feature.clinicalSignificances:`, feature.clinicalSignificances);
                    console.log(`[DEBUG]   feature.association:`, feature.association); // New debug log
                    console.log(`[DEBUG]   feature.description: "${feature.description}"`);

                    if (feature.association && feature.association.length > 0) { // Check feature.association first
                        featureDescriptionPart = ` | Association:${feature.association.map(a => a.name).join(', ')}`;
                        console.log(`[DEBUG]   Selected Association: ${featureDescriptionPart}`);
                    } else if (feature.diseaseAssociations && feature.diseaseAssociations.length > 0) {
                        featureDescriptionPart = ` | Disease:${feature.diseaseAssociations.map(da => da.diseaseName).join(', ')}`;
                        console.log(`[DEBUG]   Selected Disease Association: ${featureDescriptionPart}`);
                    } else if (feature.clinicalSignificances && feature.clinicalSignificances.length > 0) {
                        featureDescriptionPart = ` | ClinicalSignificance:${feature.clinicalSignificances.map(cs => cs.type).join(', ')}`;
                        console.log(`[DEBUG]   Selected Clinical Significance: ${featureDescriptionPart}`);
                    } else if (feature.description) {
                        featureDescriptionPart = ` | ${feature.description}`;
                        console.log(`[DEBUG]   Selected General Description: ${featureDescriptionPart}`);
                    } else {
                        console.log(`[DEBUG]   No specific description part found.`);
                    }
                    // --- END MODIFIED LOGIC ---

                    const fastaHeader = (
                        `>${modifiedCanonicalHeaderPart}` + // Use dynamic headerTag
                        ` | ${headerTag}:${mutationDescription}` + // Keep this part
                        featureDescriptionPart +
                        ` | Genomic:${firstGenomicLoc}`
                    );
                    pathogenicVariantSequences.push(`${fastaHeader}\n${modifiedSequence}`);
                }
            }
        }
        return pathogenicVariantSequences;
    }, []); // Removed currentVariantCriteria from dependencies, will pass it directly

    // Handle FASTA file upload
    const handleFileUpload = (event) => {
        const file = event.target.files[0];
        if (file) {
            setFastaFile(file);
            setResultFasta('');
            setError('');
            setStatusMessage('');
        }
    };

    // Handle variant criteria change
    const handleCriteriaChange = (event) => {
        setVariantCriteria(event.target.value);
    };

    // Process the FASTA file and generate variants
    const processFastaFile = useCallback(async () => {
        if (!fastaFile) {
            setError("Please upload a FASTA file first.");
            return;
        }

        setIsLoading(true);
        setStatusMessage("Step 1/4: Reading FASTA file...");
        setError('');
        setResultFasta('');

        try {
            const reader = new FileReader();
            reader.onload = async (e) => {
                const fastaContent = e.target.result;
                const { uniprotIds, sequences, canonicalHeaders } = getSequencesAndIdsFromFasta(fastaContent);

                if (uniprotIds.length === 0) {
                    setError("No valid UniProt IDs found in the FASTA file.");
                    setIsLoading(false);
                    return;
                }

                setStatusMessage(`Step 2/4: Downloading variant data for ${uniprotIds.length} proteins...`);
                const allPathogenicVariants = [];

                for (let i = 0; i < uniprotIds.length; i++) {
                    const uniprotId = uniprotIds[i];
                    setStatusMessage(`Step 2/4: Downloading data for ${uniprotId} (${i + 1}/${uniprotIds.length})...`);
                    const canonicalSeq = sequences[uniprotId];
                    const canonicalHeader = canonicalHeaders[uniprotId];

                    if (!canonicalSeq || !canonicalHeader) {
                        console.warn(`Warning: Canonical sequence or header not found for ${uniprotId}.`);
                        continue;
                    }

                    // ====================================================================
                    // >>>>>> KEY CHANGE HERE: REAL UNIPROT API CALL <<<<<<
                    // ====================================================================
                    let variantsData = null;
                    try {
                        const response = await fetch(`https://www.ebi.ac.uk/proteins/api/variation/${uniprotId}?format=json`);
                        if (!response.ok) {
                            // Handle HTTP errors (e.g., 404 Not Found, 500 Server Error)
                            if (response.status === 404) {
                                console.warn(`Warning: No variant data found for ${uniprotId} (404 Not Found).`);
                            } else {
                                console.error(`Error fetching data for ${uniprotId}: ${response.status} ${response.statusText}`);
                            }
                            variantsData = { features: [] }; // Assume no variants if there's an error
                        } else {
                            variantsData = await response.json();
                        }
                    } catch (fetchError) {
                        console.error(`Network or CORS error fetching data for ${uniprotId}:`, fetchError);
                        variantsData = { features: [] }; // Assume no variants if there's a network error
                    }
                    // ====================================================================
                    // END OF CHANGE

                    if (variantsData && variantsData.features && variantsData.features.length > 0) {
                        setStatusMessage(`Step 3/4: Processing variants for ${uniprotId}...`);
                        // Pass current variantCriteria to the processing function
                        const proteinPathogenicVariants = processProteinVariants(uniprotId, canonicalSeq, canonicalHeader, variantsData, variantCriteria);
                        allPathogenicVariants.push(...proteinPathogenicVariants);
                    } else {
                        // --- DEBUG LOG ---
                        console.log(`[DEBUG] No features found for UniProt ID: ${uniprotId}`);
                        // --- END DEBUG ---
                    }
                }

                setStatusMessage("Step 4/4: Combining all variants and cleaning asterisks...");
                let combinedFastaContent = allPathogenicVariants.map(entry => {
                    const [header, sequence] = entry.split('\n');
                    const cleanedSequence = sequence.replace(/\*/g, ''); // Remove asterisks
                    return `${header}\n${cleanedSequence}`;
                }).join('\n');

                if (allPathogenicVariants.length === 0) {
                    let noResultsMessage = "";
                    if (variantCriteria === 'pathogenicOnly') {
                        noResultsMessage = "No pathogenic or likely pathogenic variants were found for any of the proteins in the provided FASTA file.";
                    } else if (variantCriteria === 'allWithDescription') {
                        noResultsMessage = "No variants with descriptions were found for any of the proteins in the provided FASTA file, based on the selected criteria.";
                    } else if (variantCriteria === 'allVariants') { // New message for new option
                        noResultsMessage = "No variants (of type VARIANT or Natural variant) were found for any of the proteins in the provided FASTA file.";
                    } else {
                        noResultsMessage = "No variants were found for any of the proteins in the provided FASTA file."; // Fallback, should not happen with current options
                    }
                    setResultFasta(noResultsMessage);
                    setStatusMessage("Process completed.");
                } else {
                    setResultFasta(combinedFastaContent);
                    setStatusMessage("Process completed! Variants FASTA file generated."); // Changed to 'Variants'
                }
                setIsLoading(false);
            };
            reader.onerror = (e) => {
                setError("Error reading file: " + e.target.error);
                setIsLoading(false);
            };
            reader.readAsText(fastaFile);

        } catch (err) {
            setError("An unexpected error occurred during processing: " + err.message);
            console.error(err);
            setIsLoading(false);
        }
    }, [fastaFile, getSequencesAndIdsFromFasta, processProteinVariants, variantCriteria]); // Added variantCriteria to dependencies

    // Function to download the resulting FASTA file
    const handleDownload = () => {
        const blob = new Blob([resultFasta], { type: 'text/fasta' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'all_pathogenic_variants_combined.fasta';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    };

    // User Guide content (from the provided Markdown, converted to plain text for download)
    const userGuideContent = `# User Guide: Pathogenic Variant Generator

Welcome to the Pathogenic Variant Generator! This powerful tool helps you identify and automatically obtain protein variant sequences that are highly relevant for your research, especially those classified as pathogenic or those with detailed descriptions.

## What does this application do?

Imagine you have a list of proteins (in a FASTA file) and you need to quickly determine if they harbor genetic variants that might be associated with diseases or are otherwise of significant interest. This application streamlines that process by performing the following key functions:

1. Reads your protein file: You simply provide your FASTA file, which contains the "original" or canonical sequences of your proteins.

2. Automatically searches for variants: The application intelligently connects to a vast protein database (specifically, UniProt) and efficiently searches for all known variants associated with each of your provided proteins.

3. Filters for what's important: You have full control! The tool allows you to precisely choose which type of variants you are most interested in:

   * Pathogenic/Likely Pathogenic Variants Only: This is the strictest and most focused option. It will exclusively include variants that UniProt has explicitly classified as "Pathogenic" or "Likely pathogenic" based on their clinical significance.

   * All Variants with Descriptions (if available): This option provides a broader view. It will include any variant that comes with a textual description from UniProt, regardless of its clinical significance. This is incredibly useful for exploring variants with potential functional information.

   * All Variants (regardless of description or pathogenicity): This is the most inclusive option. It will generate sequences for every single variant feature that UniProt has registered for your proteins, without applying any additional filtering based on clinical impact or descriptive text.

4. Creates new sequences: For each variant that successfully meets your selected criteria, the application will "modify" the original, canonical sequence of your protein to accurately create the corresponding variant sequence.

5. Provides you with a new file: Finally, all the identified and generated variant sequences are neatly combined into a single, comprehensive FASTA file. This file is then ready for you to download and seamlessly integrate into other bioinformatics programs or analyses.

## How to use the application?

Using the Pathogenic Variant Generator is remarkably simple and intuitive. Just follow these straightforward steps:

### Step 1: Upload your FASTA file

1. Locate the designated area labeled "1. Upload your FASTA file (.fasta, .txt)".

2. Click on the "Choose File" (or similar) button to open your file explorer.

3. Navigate through your computer's directories and select the FASTA file that contains the protein sequences you wish to analyze.

   * Important: Please ensure that your file has either a .fasta or .txt extension.

   * Once selected, the name of your chosen file will be displayed directly below the upload button, confirming your selection.

### Step 2: Select variant criteria

1. Move to the section titled "2. Select Variant Inclusion Criteria". Here, you will define the specific types of variants you want the application to search for and generate.

2. Carefully choose one of the following radio button options:

   * Pathogenic/Likely Pathogenic Variants Only: Select this if your primary interest lies in variants with confirmed or highly probable clinical significance related to disease.

   * All Variants with Descriptions (if available): Choose this to get a wider range of variants, focusing on those for which UniProt provides additional descriptive context.

   * All Variants (regardless of description or pathogenicity): Opt for this if you require an exhaustive list of all known variants, without any filtering based on clinical impact or descriptive text.

### Step 3: Generate the variants

1. Once your file is uploaded and your variant criteria are selected, proceed by clicking the prominent "Generate Variants" button.

2. The application will then display a "Processing..." message, accompanied by a status bar (e.g., "Step 1/4: Reading FASTA file..."). This indicates the progress of the operation. Please note that the processing time may vary depending on the number of proteins in your input file and the speed of your internet connection, as the application actively queries the UniProt database.

### Step 4: Review and download the results

1. Upon successful completion of the processing, a new section labeled "Results" will appear.

2. If variants were found and generated according to your specified criteria, you will see a large text box containing the entire content of the combined FASTA file of all the new variant sequences.

3. To save this valuable data, simply click the "Download Combined FASTA" button. The file will be saved to your computer, ready for immediate use in other bioinformatics tools or for further analysis.

4. In the event that no variants are found based on your selected criteria, the application will clearly indicate this with an informative message.

That's it! With these simple steps, you will be able to efficiently use the Pathogenic Variant Generator for your research. Should you have any questions or encounter any problems during its use, please feel free to consult the technical documentation or reach out for assistance`;

    // Function to download the user guide
    const handleDownloadUserGuide = () => {
        const blob = new Blob([userGuideContent], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'Pathogenic_Variant_Generator_User_Guide.txt';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    };

    const creationDate = "July 10, 2025"; // Hardcoded creation date

    return (
        <div
            className="min-h-screen flex flex-col items-center justify-center p-4 font-sans text-gray-800"
            style={{
                background: 'linear-gradient(135deg, #f0f4f8 0%, #d9e2ec 100%)' // Subtle gradient background
            }}
        >
            <div className="relative bg-white p-8 rounded-2xl shadow-2xl w-full max-w-2xl
                            transform transition-all duration-300 hover:scale-[1.01] hover:shadow-3xl
                            border border-gray-200 z-10">
                <div className="flex flex-col items-center mb-6">
                    {/* Logo SVG */}
                    <div className="w-24 h-24 mb-4"> {/* Contenedor para el SVG */}
                        <svg width="100%" height="100%" viewBox="0 0 100 100" fill="none" xmlns="http://www.w3.org/2000/svg">
                          <rect width="100" height="100" rx="20" fill="#E0E7FF"/>
                          <text x="50" y="50" fontFamily="Inter, sans-serif" fontSize="38" fontWeight="bold" fill="#4F46E5" textAnchor="middle" alignmentBaseline="middle">PVG</text>
                          <path d="M25 65 L50 75 L75 65" stroke="#10B981" strokeWidth="4" strokeLinecap="round" fill="none"/>
                          <path d="M75 65 L70 60 M75 65 L70 70" stroke="#10B981" strokeWidth="4" strokeLinecap="round" strokeLinejoin="round"/>
                          <text x="50" y="88" fontFamily="Inter, sans-serif" font-size="12" fill="#6B7280" text-anchor="middle" alignment-baseline="middle">GENERATOR</text>
                        </svg>
                    </div>
                    <h1 className="text-4xl md:text-5xl font-extrabold text-center text-gray-900 mb-2 tracking-tight">
                        Pathogenic Variant Generator
                    </h1>
                    <p className="text-gray-600 text-center text-lg md:text-xl">
                        Explore and generate pathogenic protein variants.
                    </p>
                    {/* Botón para descargar la guía de usuario */}
                    <button
                        onClick={handleDownloadUserGuide}
                        className="mt-4 py-2 px-4 rounded-xl bg-gray-200 text-gray-700 font-semibold text-sm
                                    hover:bg-gray-300 active:bg-gray-400 shadow-sm transition-all duration-200"
                    >
                        Download User Guide
                    </button>
                </div>
                
                <p className="text-gray-700 text-center mb-8 text-base leading-relaxed">
                    Upload your canonical FASTA sequence file to identify and generate the sequences of their pathogenic variants.
                    Our tool streamlines the process, providing accurate results for your research.
                </p>

                {/* New section for Example of Utility - Now collapsible */}
                <div className="mb-6 border-b border-gray-200 py-6 px-6 rounded-lg bg-blue-50 shadow-inner">
                    <div className="flex justify-between items-center cursor-pointer" onClick={() => setShowUtilityExample(!showUtilityExample)}>
                        <h2 className="text-lg font-semibold text-blue-800">Example of Utility</h2>
                        <button className="text-blue-600 hover:text-blue-800 transition-colors">
                            {showUtilityExample ? (
                                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M5 15l7-7 7 7"></path></svg>
                            ) : (
                                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M19 9l-7 7-7-7"></path></svg>
                            )}
                        </button>
                    </div>
                    {showUtilityExample && (
                        <div className="mt-3">
                            <p className="text-sm text-blue-700 leading-relaxed">
                                Imagine you are studying a set of proteins involved in a specific disease, such as **Amyloidosis**. You have their canonical FASTA sequences.
                                Using this tool, you can upload your FASTA file, select "Pathogenic/Likely Pathogenic Variants Only", and the application will:
                            </p>
                            <ul className="list-disc list-inside mt-2 space-y-1 text-sm text-blue-700 leading-relaxed">
                                <li>Automatically query UniProt for each protein.</li>
                                <li>Identify variants like **TTR p.L48M** or **APOE p.R130C**, which are known to be pathogenic for Amyloidosis or Alzheimer's disease, respectively.</li>
                                <li>Generate new FASTA entries for these specific variants, **including pathology annotation in their headers**, allowing you to easily use them for downstream analysis, such as structural modeling, functional studies, **proteomic studies (e.g., mass spectrometry-based proteomics)**, or **integrating into search engines for variant identification**, to understand their impact.</li>
                            </ul>
                            <p className="text-sm text-blue-700 leading-relaxed mt-2">
                                This streamlines the process of obtaining mutated protein sequences that are directly relevant to your research on disease mechanisms.
                            </p>
                        </div>
                    )}
                </div>

                <div className="mb-6 border-t border-b border-gray-200 py-6 px-6 rounded-lg bg-gray-50 shadow-inner"> {/* Added px-6 */}
                    <label htmlFor="fasta-upload" className="block text-gray-800 text-lg font-semibold mb-3">
                        1. Upload your FASTA file (.fasta, .txt)
                    </label>
                    <input
                        type="file"
                        id="fasta-upload"
                        accept=".fasta,.txt"
                        onChange={handleFileUpload}
                        className="block w-full text-base text-gray-700
                                file:mr-4 file:py-2.5 file:px-5
                                file:rounded-full file:border-0
                                file:text-base file:font-semibold
                                file:bg-blue-100 file:text-blue-700
                                hover:file:bg-blue-200 cursor-pointer rounded-lg border border-gray-300 p-2.5 transition-colors duration-200"
                    />
                    {fastaFile && (
                        <p className="mt-3 text-sm text-gray-600">Selected file: <span className="font-medium text-blue-700">{fastaFile.name}</span></p>
                    )}
                </div>

                {/* New section for variant criteria selection */}
                <div className="mb-6 border-b border-gray-200 py-6 px-6 rounded-lg bg-gray-50 shadow-inner mt-6"> {/* Added px-6 and mt-6 */}
                    <label className="block text-gray-800 text-lg font-semibold mb-3">
                        2. Select Variant Inclusion Criteria
                    </label>
                    <div className="flex flex-col space-y-3">
                        <div className="flex items-center">
                            <input
                                type="radio"
                                id="pathogenicOnly"
                                name="variantCriteria"
                                value="pathogenicOnly"
                                checked={variantCriteria === 'pathogenicOnly'}
                                onChange={handleCriteriaChange}
                                className="form-radio h-5 w-5 text-blue-600 border-gray-300 focus:ring-blue-500"
                            />
                            <label htmlFor="pathogenicOnly" className="ml-3 text-gray-700 text-base">
                                Pathogenic/Likely Pathogenic Variants Only
                            </label>
                        </div>
                        {/* Explanation for Pathogenic/Likely Pathogenic Variants Only */}
                        <p className="ml-8 text-sm text-gray-500 -mt-2">
                            Includes only variants explicitly classified as "Pathogenic" or "Likely pathogenic" by UniProt, based on clinical significance.
                        </p>

                        <div className="flex items-center">
                            <input
                                type="radio"
                                id="allWithDescription"
                                name="variantCriteria"
                                value="allWithDescription"
                                checked={variantCriteria === 'allWithDescription'}
                                onChange={handleCriteriaChange}
                                className="form-radio h-5 w-5 text-blue-600 border-gray-300 focus:ring-blue-500"
                            />
                            <label htmlFor="allWithDescription" className="ml-3 text-gray-700 text-base">
                                All Variants with Descriptions (if available)
                            </label>
                        </div>
                        {/* Explanation for All Variants with Descriptions */}
                        <p className="ml-8 text-sm text-gray-500 -mt-2">
                            Includes all variants (both "VARIANT" and "Natural variant" types) that have a descriptive text provided by UniProt, regardless of their clinical significance.
                        </p>

                        <div className="flex items-center">
                            <input
                                type="radio"
                                id="allVariants"
                                name="variantCriteria"
                                value="allVariants"
                                checked={variantCriteria === 'allVariants'}
                                onChange={handleCriteriaChange}
                                className="form-radio h-5 w-5 text-blue-600 border-gray-300 focus:ring-blue-500"
                            />
                            <label htmlFor="allVariants" className="ml-3 text-gray-700 text-base">
                                All Variants (regardless of description or pathogenicity)
                            </label>
                        </div>
                        {/* Explanation for All Variants */}
                        <p className="ml-8 text-sm text-gray-500 -mt-2">
                            Includes every variant feature (both "VARIANT" and "Natural variant" types) returned by UniProt for the protein, without any filtering based on description or clinical significance.
                        </p>
                    </div>
                </div>

                <div className="mb-6 text-center">
                    <button
                        onClick={processFastaFile}
                        disabled={!fastaFile || isLoading}
                        className={`w-full py-3.5 px-6 rounded-xl text-white font-bold text-xl
                                    transition-all duration-300 ease-in-out
                                    ${!fastaFile || isLoading
                                        ? 'bg-gray-400 cursor-not-allowed opacity-75'
                                        : 'bg-blue-800 hover:bg-blue-900 active:bg-blue-950 shadow-lg hover:shadow-xl' // Darker blue button style
                                    }`}
                    >
                        {isLoading ? (
                            <span className="flex items-center justify-center">
                                <svg className="animate-spin -ml-1 mr-3 h-6 w-6 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                                </svg>
                                {statusMessage || 'Processing...'}
                            </span>
                        ) : (
                            'Generate Variants' // Changed button text to be more general
                        )}
                    </button>
                    {statusMessage && !isLoading && (
                        <p className="mt-4 text-blue-700 font-medium text-base">{statusMessage}</p>
                    )}
                    {error && (
                        <p className="mt-4 text-red-600 font-medium text-base bg-red-50 p-3 rounded-lg border border-red-200">{error}</p>
                    )}
                </div>

                {resultFasta && (
                    <div className="mt-8 pt-6 border-t border-gray-200">
                        <h2 className="text-2xl font-bold text-gray-800 mb-4 text-center">Results</h2>
                        <div className="bg-gray-50 p-4 rounded-lg border border-gray-200 mb-4 max-h-60 overflow-y-auto shadow-inner">
                            <pre className="text-sm text-gray-800 whitespace-pre-wrap break-words">{resultFasta}</pre>
                        </div>
                        <div className="text-center">
                            <button
                                onClick={handleDownload}
                                className="py-2.5 px-6 rounded-xl bg-green-600 text-white font-bold text-lg
                                            hover:bg-green-700 active:bg-green-800 shadow-md" // Flatter button style
                            >
                                Download Combined FASTA
                            </button>
                        </div>
                    </div>
                )}
            </div>
            {/* Footer con información de afiliación y fecha de creación */}
            <p className="mt-8 text-gray-500 text-xs text-center max-w-xl leading-relaxed z-10">
                Note: This version attempts to make direct calls to the UniProt API. If you experience errors, it might be due to CORS restrictions in your environment. For production environments, consider a backend proxy.
                <br />
                <span className="font-semibold text-gray-600 mt-2 block">
                    Dr. Carlos Fuentes, Proteomics Unit-SCAI, UCO
                    <br />
                    Email: <a href="mailto:proteomica2@uco.es" className="text-blue-500 hover:underline">proteomica2@uco.es</a>
                    <br />
                    ORCID: <a href="https://orcid.org/0000-0001-9940-9065" target="_blank" rel="noopener noreferrer" className="text-blue-500 hover:underline">0000-0001-9940-9065</a>
                </span>
                <br />
                <span className="text-gray-500 mt-2 block">
                    Created: {creationDate}
                </span>
            </p>
        </div>
    );
}

export default App;
