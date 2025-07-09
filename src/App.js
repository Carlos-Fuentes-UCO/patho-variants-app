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

                const match = trimmedLine.match(/^>[a-z]{2}\|([A-Z0-9]+)\|/);
                if (match) {
                    currentId = match[1];
                } else {
                    const parts = trimmedLine.substring(1).split(' ');
                    if (parts.length > 0) {
                        const idMatchFallback = parts[0].match(/\|([A-Z0-9]+)\|/);
                        if (idMatchFallback) {
                            currentId = idMatchFallback[1];
                        } else {
                            currentId = parts[0];
                        }
                    } else {
                        currentId = `UNKNOWN_ID_${Math.random().toString(36).substring(7)}`; // Generate random ID
                    }
                }

                if (currentId && !uniprotIds.includes(currentId)) {
                    uniprotIds.push(currentId);
                }
                canonicalHeaders[currentId] = trimmedLine; // Save the complete header line
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
                if ((feature.type === "VARIANT" || feature.type === "Natural variant") && feature.description) {
                    includeFeature = true;
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Included variant for 'All Variants with Descriptions': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Description: "${feature.description}"`);
                    // --- END DEBUG LOG ---
                } else {
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Skipped variant for 'All Variants with Descriptions': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Has Description: ${!!feature.description}`);
                    // --- END DEBUG LOG ---
                }
            } else if (currentVariantCriteria === 'allVariants') {
                // Include all 'VARIANT' or 'Natural variant' features
                if (feature.type === "VARIANT" || feature.type === "Natural variant") {
                    includeFeature = true;
                    // --- DEBUG LOG ---
                    console.log(`[DEBUG] Included variant for 'All Variants': UniProt ID: ${uniprotId}, Feature Type: ${feature.type}, Has Description: ${!!feature.description}`);
                    // --- END DEBUG LOG ---
                }
            }

            if (includeFeature) {
                const beginPos = parseInt(feature.begin);
                const endPos = parseInt(feature.end);
                
                // Be more robust with wildType/mutatedType, falling back to alternativeSequence if direct fields are missing
                let wildType = feature.wildType || (feature.alternativeSequence && feature.alternativeSequence.originalSequence);
                let mutatedType = feature.mutatedType || (feature.alternativeSequence && feature.alternativeSequence.alternativeSequences && feature.alternativeSequence.alternativeSequences[0]);

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
                    const currentAminoAcidInFasta = variantSequenceList[beginPos - 1];

                    if (wildType && mutatedType) { // Substitution
                        if (currentAminoAcidInFasta !== wildType) {
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
                    let modifiedCanonicalHeaderPart = cleanedOriginalHeader;

                    // Find the second '|' to insert ftId (e.g., after P06396 in ">sp|P06396|GELS_HUMAN...")
                    const firstPipeIndex = cleanedOriginalHeader.indexOf('|');
                    const secondPipeIndex = cleanedOriginalHeader.indexOf('|', firstPipeIndex + 1);

                    if (secondPipeIndex !== -1 && feature.ftId) {
                        // Reconstruct the header with ftId
                        modifiedCanonicalHeaderPart = 
                            cleanedOriginalHeader.substring(0, secondPipeIndex) +
                            '-' + feature.ftId +
                            cleanedOriginalHeader.substring(secondPipeIndex);
                    }

                    // Add feature.description if it exists
                    const featureDescriptionPart = feature.description ? ` | ${feature.description}` : '';

                    const fastaHeader = (
                        `>${modifiedCanonicalHeaderPart} ` +
                        `| PATHOGENIC_VARIANT:${mutationDescription}` +
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
                        // --- END DEBUG LOG ---
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

    return (
        <div
            className="min-h-screen flex flex-col items-center justify-center p-4 font-sans text-gray-800"
            style={{
                backgroundColor: '#ffffff' // Pure white background
            }}
        >
            <div className="relative bg-white p-8 rounded-2xl shadow-lg w-full max-w-2xl
                            border border-gray-200 z-10"> {/* Removed hover effects */}
                <div className="flex flex-col items-center mb-6">
                    {/* No Logo - Title is prominent */}
                    <h1 className="text-4xl md:text-5xl font-extrabold text-center text-gray-900 mb-2 tracking-tight">
                        Pathogenic Variant Generator
                    </h1>
                    <p className="text-gray-600 text-center text-lg md:text-xl">
                        Explore and generate pathogenic protein variants.
                    </p>
                </div>
                
                <p className="text-gray-700 text-center mb-8 text-base leading-relaxed">
                    Upload your canonical FASTA sequence file to identify and generate the sequences of their pathogenic variants.
                    Our tool streamlines the process, providing accurate results for your research.
                </p>

                <div className="mb-6 border-t border-b border-gray-200 py-6 px-4 rounded-lg bg-gray-50 shadow-inner">
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
                <div className="mb-6 border-b border-gray-200 py-6 px-4 rounded-lg bg-gray-50 shadow-inner">
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
                                        : 'bg-blue-600 hover:bg-blue-700 active:bg-blue-800 shadow-md' // Flatter button style
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
            <p className="mt-8 text-gray-500 text-xs text-center max-w-xl leading-relaxed z-10">
                Note: This version attempts to make direct calls to the UniProt API. If you experience errors, it might be due to CORS restrictions in your environment. For production environments, consider a backend proxy.
            </p>
        </div>
    );
}

export default App;
