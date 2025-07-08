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
    const processProteinVariants = useCallback((uniprotId, canonicalSequence, canonicalHeaderLine, variantsData) => {
        const pathogenicVariantSequences = [];
        if (!variantsData || !variantsData.features) {
            return [];
        }

        for (const feature of variantsData.features) {
            if (feature.type === "VARIANT") {
                const beginPos = parseInt(feature.begin);
                const endPos = parseInt(feature.end);
                const wildTypeAa = feature.wildType;
                const mutatedTypeAa = feature.mutatedType;

                let mutationDescription = "";
                if (wildTypeAa && mutatedTypeAa) {
                    mutationDescription = `p.${wildTypeAa}${beginPos}${mutatedTypeAa}`;
                } else if (wildTypeAa && !mutatedTypeAa && beginPos === endPos) {
                    mutationDescription = `p.del${wildTypeAa}${beginPos}`;
                } else if (wildTypeAa && !mutatedTypeAa && beginPos !== endPos) {
                    mutationDescription = `p.del${wildTypeAa}${beginPos}-${endPos}`;
                } else if (!wildTypeAa && mutatedTypeAa) {
                    mutationDescription = `p.ins{mutatedTypeAa}${beginPos}`;
                }

                if (feature.description) {
                    mutationDescription += ` (${feature.description})`;
                } else if (!mutationDescription) {
                    mutationDescription = `Variant at pos ${beginPos}`;
                }

                let isPathogenic = false;
                if (feature.clinicalSignificances) {
                    for (const cs of feature.clinicalSignificances) {
                        if (cs.type === "Pathogenic" || cs.type === "Likely pathogenic") {
                            isPathogenic = true;
                            break;
                        }
                    }
                }

                // Declare modifiedSequence here, outside the isPathogenic block
                // This ensures it's always declared for each feature iteration.
                let modifiedSequence = ''; 

                if (isPathogenic) {
                    let variantSequenceList = Array.from(canonicalSequence);
                    
                    // Apply the mutation
                    let mutationAppliedSuccessfully = false;

                    if (beginPos - 1 < 0 || beginPos - 1 >= variantSequenceList.length) {
                        // Position out of bounds, mutation not applied
                    } else if (wildTypeAa && mutatedTypeAa) { // Substitution
                        if (variantSequenceList[beginPos - 1] === wildTypeAa) {
                            variantSequenceList[beginPos - 1] = mutatedTypeAa;
                            mutationAppliedSuccessfully = true;
                        }
                    } else if (wildTypeAa && !mutatedTypeAa) { // Deletion
                        if (endPos <= variantSequenceList.length && beginPos <= endPos) {
                            const expectedDeleted = canonicalSequence.substring(beginPos - 1, endPos);
                            if (expectedDeleted === wildTypeAa) {
                                variantSequenceList.splice(beginPos - 1, endPos - (beginPos - 1));
                                mutationAppliedSuccessfully = true;
                            }
                        }
                    } else if (!wildTypeAa && mutatedTypeAa) { // Insertion
                        if (beginPos - 1 <= variantSequenceList.length) {
                            variantSequenceList.splice(beginPos - 1, 0, mutatedTypeAa);
                            mutationAppliedSuccessfully = true;
                        }
                    }

                    if (mutationAppliedSuccessfully) {
                        modifiedSequence = variantSequenceList.join(''); // Assign only if mutation was applied
                        
                        const genomicLocInfo = feature.genomicLocation || [''];
                        const firstGenomicLoc = genomicLocInfo[0] || '';

                        // Use the standalone lstripString function
                        const cleanedCanonicalHeader = lstripString(canonicalHeaderLine.trim(), '>');
                        const fastaHeader = (
                            `>${cleanedCanonicalHeader} ` +
                            `| PATHOGENIC_VARIANT:${mutationDescription} ` +
                            `| Genomic:${firstGenomicLoc}`
                        );
                        pathogenicVariantSequences.push(`${fastaHeader}\n${modifiedSequence}`);
                    }
                }
            }
        }
        return pathogenicVariantSequences;
    }, []);

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
                        const proteinPathogenicVariants = processProteinVariants(uniprotId, canonicalSeq, canonicalHeader, variantsData);
                        allPathogenicVariants.push(...proteinPathogenicVariants);
                    } else {
                        // console.log(`No variants with 'features' found for ${uniprotId}.`);
                    }
                }

                setStatusMessage("Step 4/4: Combining all variants and cleaning asterisks...");
                let combinedFastaContent = allPathogenicVariants.map(entry => {
                    const [header, sequence] = entry.split('\n');
                    const cleanedSequence = sequence.replace(/\*/g, ''); // Remove asterisks
                    return `${header}\n${cleanedSequence}`;
                }).join('\n');

                if (allPathogenicVariants.length === 0) {
                    setResultFasta("No pathogenic variants were found for any of the proteins in the provided FASTA file.");
                    setStatusMessage("Process completed.");
                } else {
                    setResultFasta(combinedFastaContent);
                    setStatusMessage("Process completed! Pathogenic variants FASTA file generated.");
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
    }, [fastaFile, getSequencesAndIdsFromFasta, processProteinVariants]);

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
                background: 'linear-gradient(135deg, #f0f4f8 0%, #d9e2ec 100%)' // Subtle gradient background
            }}
        >
            <div className="relative bg-white p-8 rounded-2xl shadow-2xl w-full max-w-2xl
                            transform transition-all duration-300 hover:scale-[1.01] hover:shadow-3xl
                            border border-gray-200 z-10">
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

                <div className="mb-6 text-center">
                    <button
                        onClick={processFastaFile}
                        disabled={!fastaFile || isLoading}
                        className={`w-full py-3.5 px-6 rounded-xl text-white font-bold text-xl
                                    transition-all duration-300 ease-in-out transform
                                    ${!fastaFile || isLoading
                                        ? 'bg-gray-400 cursor-not-allowed opacity-75'
                                        : 'bg-gradient-to-r from-blue-700 to-indigo-800 hover:from-blue-800 hover:to-indigo-900 active:from-blue-900 active:to-indigo-950 shadow-lg hover:shadow-xl'
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
                            'Generate Pathogenic Variants'
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
                                className="py-2.5 px-6 rounded-xl bg-gradient-to-r from-green-600 to-teal-700 text-white font-bold text-lg
                                            hover:from-green-700 hover:to-teal-800 active:from-green-800 active:to-teal-900 shadow-md hover:shadow-lg
                                            transition-all duration-300 ease-in-out transform"
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