import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

pipeline {
    agent any
    // Define parameters that can be edited via the Jenkins UI
    parameters {
        booleanParam(name: 'TRIGGER_BUILD', defaultValue: true, description: 'Check this to trigger the build. If unchecked, the build will not be triggered and only the config.json will be generated.')
        booleanParam(name: 'CREATE_CELLDB_METADATA', defaultValue: true, description: 'Check this to trigger the create_celldb_metadata job.')
        booleanParam(name: 'PULL_POOL_ID', defaultValue: false, description: 'Flag indicating whether to pull pool IDs from CellDB - only applicable to cell sets (i.e. EXT.PR500.CS01.1.A, EXT.PR500.CS01.1.B, etc).')
        booleanParam(name: 'COLLATE_FASTQ_READS', defaultValue: true, description: 'Check this to trigger the collate_fastq_reads job.')
        booleanParam(name: 'FILTER_COUNTS', defaultValue: true, description: 'Check this to trigger the filter_counts job.')
        booleanParam(name: 'REMOVE_DATA', defaultValue: false, description: 'Select if there is experimental data that needs to be removed before normalization. TODO: expand on this.')
        booleanParam(name: 'CBNORMALIZE', defaultValue: true, description: 'Check this to trigger the CBnormalize job.')
        booleanParam(name: 'COMPUTE_LFC', defaultValue: true, description: 'Check this to trigger the compute_l2fc job.')
        booleanParam(name: 'COLLAPSE', defaultValue: true, description: 'Check this to trigger the collapse job.')
        booleanParam(name: 'QC_IMAGES', defaultValue: true, description: 'Check this to trigger the QC job.')
        booleanParam(name: 'JOIN_METADATA', defaultValue: true, description: 'Check this to trigger the join_metadata job. This should be checked you are planning on ruunning CONVERT_SUSHI.')
        booleanParam(name: 'CONVERT_SUSHI', defaultValue: false, description: 'Convert output column headers to format for MTS pipeline and upload to s3.')
        booleanParam(name: 'RUN_EPS_QC', defaultValue: false, description: 'Run EPS QC')
        string(name: 'BUILD_DIR', defaultValue: '/cmap/obelix/pod/prismSeq/', description: 'Output path to deposit build. Format should be /directory/PROJECT_CODE/BUILD_NAME')
        string(name: 'BUILD_NAME', defaultValue: '', description: 'Build name')
        string(name: 'SCREEN', defaultValue: '', description: 'Screen name from COMET, necessary if using COMET for sample metadata.')
        string(name: 'SEQ_TYPE', defaultValue: 'DRAGEN', description: 'Choose DRAGEN, MiSeq, HiSeq, or NovaSeq. MiSeq and HiSeq/NovaSeq return files named differently. This setting sets the INDEX_1, INDEX_2, and BARCODE_SUFFIX parameters in fastq2readcount. Select DRAGEN if fastq files are from the DRAGEN pipeline from GP. Choosing NovaSeq reverses index 2.')
        string(name: 'DAYS', defaultValue: '', description: 'If running the sushi_to_mts module, provide any days/timepoints (separated by commas) that should be dropped from output data. No quotes needed (ie, 2,8).')

        // pipeline version
        string(name: 'GIT_BRANCH', defaultValue: 'main', description: 'Pipeline branch to use')
        booleanParam(name: 'USE_LATEST', defaultValue: true, description: 'Check this to use the most up to date version from the specified branch. If not checked, will use the specified commit.')
        string(name: 'COMMIT_ID', defaultValue: '', description: 'Specific commit ID to use (leave empty if using the latest commit in the branch or if already specified in the config file.)')

        // Most common parameters
        string(name: 'SEQUENCING_INDEX_COLS', defaultValue: 'flowcell_names,index_1,index_2', description: 'Used in COLLATE_FASTQ_READS, this is a comma separate list of sequencing columns in the sample meta that are need to identify every PCR well in the run.')
        string(name: 'ID_COLS', defaultValue: 'pcr_plate,pcr_well', description: 'Used in COLLATE_FASTQ_READS, columns to concat to create unique ID for each sample-replicate')
        string(name: 'CELL_LINE_COLS', defaultValue: 'DepMap_ID', description: 'Columns in intermediate files that describe a read or cell line')
        string(name: 'SIG_COLS', defaultValue: 'cell_set,treatment,dose,dose_unit,day', description: 'Signature columns, these describe unique treatment conditions and generally should not include replicate information.')
        string(name: 'CONTROL_COLS', defaultValue: 'cell_set,day', description: 'Set of columns that define individual negative control conditions.')

        // Sushi Input files
        string(name: 'RAW_COUNTS_UNCOLLAPSED', defaultValue: 'raw_counts_uncollapsed.csv', description: 'Filename in BUILD_DIR containing nori output')
        string(name: 'SAMPLE_META', defaultValue: 'sample_meta.csv', description: 'File name of sample metadata within the BUILD_DIR directory.')
        string(name: 'CELL_SET_META', defaultValue: 'cell_set_meta.csv', description: 'Cell Set Metadata. Static cell_line_meta location: /data/vdb/prismSeq/cell_set_meta.csv')
        string(name: 'CELL_LINE_META', defaultValue: 'cell_line_meta.csv', description: 'File in BUILD_DIR containing cell line metadata')
        string(name: 'CONTROL_BARCODE_META', defaultValue: 'CB_meta.csv', description: 'Metadata for control barcodes.')
        string(name: 'ASSAY_POOL_META', defaultValue: 'assay_pool_meta.txt', description: 'File in BUILD_DIR containing assay pool metadata')

        // Additional parameters ordered by when they first appear
        string(name: 'BARCODE_COL', defaultValue: 'forward_read_cl_barcode', description: 'Used in COLLATE_FASTQ_READS, the name of the column containing the read')
        string(name: 'LOW_ABUNDANCE_THRESHOLD', defaultValue: '20', description: 'Used in COLLATE_FASTQ_READS, threshold for unknown barcodes')
        string(name: 'PSEUDOCOUNT', defaultValue: '20', description: 'Used in CBNORMALIZE, the pesudocount value for log transformations.')
        string(name: 'COUNT_COL_NAME', defaultValue: 'normalized_n', description: 'Used in COMPUTE_LFC, the name of the numeric column to use for calculations')
        string(name: 'CTL_TYPES', defaultValue: 'negcon', description: 'Used in COMPUTE_LFC, the value in trt_type that indicates the negative controls')
        string(name: 'COUNT_THRESHOLD', defaultValue: '40', description: 'Used in COMPUTE_LFC, the count threshold for the collapsed negative controls. Cell lines in the negative controls below this threshold will be dropped from log2 fold change calculations.')

        // Files created by sushi
        string(name: 'PRISM_BARCODE_COUNTS', defaultValue: 'prism_barcode_counts.csv', description: 'Filename in BUILD_DIR containing PRISM barcode counts')
        string(name: 'UNKNOWN_BARCODE_COUNTS', defaultValue: 'unknown_barcode_counts.csv', description: 'Filename in BUILD_DIR containing unknown barcode counts')
        string(name: 'ANNOTATED_COUNTS', defaultValue: 'annotated_counts.csv', description: 'Filename in BUILD_DIR containing annotated counts')
        string(name: 'FILTERED_COUNTS', defaultValue: 'filtered_counts.csv', description: 'Filename in BUILD_DIR containing filtered counts')
        string(name: 'NORMALIZED_COUNTS', defaultValue: 'normalized_counts.csv', description: 'Filename in BUILD_DIR containing normalized counts')
        string(name: 'LFC', defaultValue: 'l2fc.csv', description: 'Filename containing log2 fold change values')
        string(name: 'COLLAPSED_LFC', defaultValue: 'collapsed_l2fc.csv', description: 'Filename in BUILD_DIR containing replicate collapsed l2fc values')
        // Other
        string(name: 'API_URL', defaultValue: 'https://api.clue.io/api/', description: 'API URL')
    }

    environment {
        CONFIG_FILE_PATH = "${params.BUILD_DIR}/config.json"
    }

    stages {
        stage('Checkout') {
            steps {
                script {
                    def config = [:]
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    if (params.USE_LATEST) {
                        // Checkout the latest commit from the specified branch
                        checkout([$class: 'GitSCM',
                                  branches: [[name: "*/${params.GIT_BRANCH}"]],
                                  doGenerateSubmoduleConfigurations: false,
                                  extensions: [],
                                  userRemoteConfigs: scm.userRemoteConfigs
                        ])
                        // Overwrite the commit ID in the config with the latest commit
                        def latestCommitID = sh(script: 'git rev-parse --short HEAD', returnStdout: true).trim()
                        config.COMMIT = latestCommitID
                        echo "Using latest commit: ${latestCommitID}"
                    } else {
                        // Use the commit ID specified in the config.json
                        if (config.COMMIT) {
                            checkout([$class: 'GitSCM',
                                      branches: [[name: config.COMMIT]],
                                      doGenerateSubmoduleConfigurations: false,
                                      extensions: [[$class: 'RelativeTargetDirectory', relativeTargetDir: '']],
                                      submoduleCfg: [],
                                      userRemoteConfigs: scm.userRemoteConfigs
                            ])
                            echo "Using commit from config.json: ${config.COMMIT}"
                        } else {
                            error("COMMIT not specified in config.json and USE_LATEST is false.")
                        }
                    }
                }
            }
        }

        stage('Generate JSON Config') {
            steps {
                script {
                    def paramList = [
                        'SEQ_TYPE', 'API_URL', 'BUILD_DIR', 'INDEX_1', 'INDEX_2', 'BARCODE_SUFFIX',
                        'BUILD_NAME', 'CONVERT_SUSHI', 'PULL_POOL_ID', 'RUN_EPS_QC', 'REMOVE_DATA', 'DAYS',
                        'COUNTS', 'SCREEN',

                        // metadata files
                        'SAMPLE_META', 'CELL_SET_META', 'CELL_LINE_META', 'CONTROL_BARCODE_META', 'ASSAY_POOL_META',

                        // sushi files
                        'RAW_COUNTS_UNCOLLAPSED', 'PRISM_BARCODE_COUNTS', 'UNKNOWN_BARCODE_COUNTS', 
                        'ANNOTATED_COUNTS', 'FILTERED_COUNTS', 'NORMALIZED_COUNTS', 'LFC', 'COLLAPSED_LFC',

                        // collate_fastq_reads parameters
                        'SEQUENCING_INDEX_COLS', 'ID_COLS', 'BARCODE_COL', 'LOW_ABUNDANCE_THRESHOLD', 'REVERSE_INDEX2',

                        // normalize parameters
                        'PSEUDOCOUNT',

                        // compute_l2fc paramters
                        'SIG_COLS', 'CONTROL_COLS', 'CELL_LINE_COLS', 'COUNT_COL_NAME', 'CTL_TYPES', 'COUNT_THRESHOLD',
                    ]

                    def config = [:]

                    // Load existing config if it exists
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    // Add parameters from Jenkins UI only if they don't exist in the config
                    paramList.each { paramName ->
                        if (!config.containsKey(paramName) && params.containsKey(paramName)) {
                            config[paramName] = params[paramName]
                        }
                    }

                    // Add or overwrite specific keys
                    if (!config.containsKey('API_KEY')) {
                        config.API_KEY = sh(script: 'cat /local/jenkins/.clue_api_key', returnStdout: true).trim()
                    }

                    // Explicit settings that are programmatically derived
                    config.REVERSE_INDEX2 = config.SEQ_TYPE == 'DRAGEN'

                    // Write the config back to file after all updates
                    writeFile file: env.CONFIG_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                    echo "Generated config.json: ${config}"
                }
            }
        }

        stage('Add Commit ID to Config') {
            steps {
                script {
                    def commitID = sh(script: 'git rev-parse --short HEAD', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg commit "$commitID" '. + {COMMIT: \$commit}' $configFilePath > ${configFilePath}.tmp && mv ${configFilePath}.tmp $configFilePath
                    """
                    echo "Added commit ID to config.json: $commitID"
                }
            }
        }

        stage('Add Timestamp to Config') {
            steps {
                script {
                    def buildtime = sh(script: 'date -u +"%Y-%m-%dT%H:%M:%SZ"', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg buildtime "$buildtime" '. + {TIMESTAMP: \$buildtime}' $configFilePath > ${configFilePath}.tmp && mv ${configFilePath}.tmp $configFilePath
                    """
                    echo "Added build timestamp to config.json: $buildtime"
                }
            }
        }

        stage('Show available podman containers') {
            steps {
                script {
                    sh '/usr/bin/podman images'
                }
            }
        }

        stage('Run Scripts in Container') {
            steps {
                script {
                    if (params.TRIGGER_BUILD) {
                        def scriptsToRun = []
                        if (params.CREATE_CELLDB_METADATA) {
                            scriptsToRun.add('create_celldb_metadata.sh')
                        }
                        if (params.SCREEN?.trim()) {
                            scriptsToRun.add('create_sample_meta.sh')
                        }
                        if (params.COLLATE_FASTQ_READS) {
                            scriptsToRun.add('collate_fastq_reads.sh')
                        }
                        if (params.FILTER_COUNTS) {
                            scriptsToRun.add('filter_counts.sh')
                        }
                        if (params.CBNORMALIZE) {
                            scriptsToRun.add('CBnormalize.sh')
                        }
                        if (params.COMPUTE_LFC) {
                            scriptsToRun.add('compute_l2fc.sh')
                        }
                        if (params.COLLAPSE) {
                            scriptsToRun.add('collapse_replicates.sh')
                        }
                        if (params.QC_IMAGES) {
                            scriptsToRun.add('filteredCounts_QC.sh')
                        }
                        if (params.JOIN_METADATA) {
                            scriptsToRun.add('join_metadata.sh')
                        }
                        if (params.RUN_EPS_QC) {
                            scriptsToRun.add('eps_qc.sh')
                        }
                        if (params.CONVERT_SUSHI) {
                            scriptsToRun.add('seq_to_mts.sh')
                        }

                        scriptsToRun.each { scriptName ->
                            echo "Running script: ${scriptName}"

                            sh """
                                chmod +x $WORKSPACE/scripts/launch_job.sh
                                $WORKSPACE/scripts/launch_job.sh $scriptName
                            """
                        }
                    }
                }
            }
        }
    }

    post {
        always {
            script {
                def buildDir = params.BUILD_DIR
                def buildName = params.BUILD_NAME ?: "default_build_name"
                def logFilePath = "${buildDir}/${buildName}_console_output.log"

                // Capture the entire console log
                def log = currentBuild.rawBuild.getLog(999999)

                // Write the log to the specified file
                writeFile file: logFilePath, text: log.join("\n")

                echo "Console output has been saved to ${logFilePath}"
            }
        }
        success {
            script {
                echo 'Build completed successfully.'
            }
        }
    }
}
