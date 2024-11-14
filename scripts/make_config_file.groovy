import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

String sectionHeaderStyleGreen = ' color: white; background: green; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; '
String sectionHeaderStyleRed = ' color: white; background: red; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; '
String separatorStyleCss = ' border: 0; border-bottom: 1px dashed #ccc; background: #999; '

pipeline {
    agent any
    // Define parameters that can be edited via the Jenkins UI
    parameters {
        separator(
          name: "Group_1",
          sectionHeader: "User Inputs",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleGreen
        )
        // Check boxes of modules to run
        booleanParam(name: 'TRIGGER_BUILD', defaultValue: true, description: 'Check this to trigger the build. If unchecked, the build will not be triggered and only the config.json will be generated.')
        booleanParam(name: 'CREATE_CELLDB_METADATA', defaultValue: true, description: 'Check this to trigger the create_celldb_metadata job.')
        booleanParam(name: 'CREATE_SAMPLE_META', defaultValue: false, description: 'Get metadata from COMET, use only if screen is registered.')
        string(name: 'SCREEN', defaultValue: '', description: 'If CREATE_SAMPLE_META is checked, provide the screen name from COMET.')
        booleanParam(name: 'COLLATE_FASTQ_READS', defaultValue: true, description: 'Check this to trigger the collate_fastq_reads job.')
        booleanParam(name: 'FILTER_COUNTS', defaultValue: true, description: 'Check this to trigger the filter_counts job.')
        booleanParam(name: 'REMOVE_DATA', defaultValue: false, description: 'Select if there is experimental data that needs to be removed prior to normalization.')
        booleanParam(name: 'CBNORMALIZE', defaultValue: true, description: 'Run normalization.')
        booleanParam(name: 'COMPUTE_LFC', defaultValue: true, description: 'Compute the fold changes.')
        booleanParam(name: 'COLLAPSE', defaultValue: true, description: 'Collapse replicates.')
        booleanParam(name: 'QC_IMAGES', defaultValue: true, description: 'Check this to trigger the QC images job.')
        booleanParam(name: 'CONVERT_SUSHI', defaultValue: false, description: 'Convert output column headers to format for MTS pipeline and upload to s3.')
        string(name: 'DAYS', defaultValue: '', description: 'If running the sushi_to_mts module, provide any days/timepoints (separated by commas) that should be dropped from output data. No quotes needed (ie, 2,8).')
        booleanParam(name: 'RUN_EPS_QC', defaultValue: false, description: 'Run EPS QC')

        // Parameters we expect users to change
        string(name: 'BUILD_DIR', defaultValue: '/cmap/obelix/pod/prismSeq/', description: 'Output path to deposit build. Format should be /directory/PROJECT_CODE/BUILD_NAME')
        string(name: 'BUILD_NAME', defaultValue: '', description: 'Build name')
        string(name: 'SEQ_TYPE', defaultValue: 'DRAGEN', description: 'Choose DRAGEN, MiSeq, HiSeq, or NovaSeq. MiSeq and HiSeq/NovaSeq return files named differently. This setting sets the INDEX_1, INDEX_2, and BARCODE_SUFFIX parameters in fastq2readcount. Select DRAGEN if fastq files are from the DRAGEN pipeline from GP. Choosing NovaSeq reverses index 2.')
        string(name: 'SIG_COLS', defaultValue: 'cell_set,pert_name,pert_dose,pert_dose_unit,day', description: 'List of signature columns found in the sample meta that describeunique treatment conditions.This defaults to \"cell_set,pert_name,pert_dose,pert_dose_unit,day\". Generally, this list should NOT include replicate information such as \"tech_rep\" or \"bio_rep\". This paramter is first used in COMPUTE_LFC.')
        string(name: 'CTL_TYPES', defaultValue: 'negcon', description: 'Value in the pert_type column of the sample meta that identifies the negative contols. This defaults to \"negcon\" and is used in COMPUTE_LFC.')
        string(name: 'CONTROL_COLS', defaultValue: 'cell_set,day', description: 'List of columns found in the sample meta that describe individual negative control conditions. This defaults to \"cell_set,day\" and can be expanded to include \"pert_vehicle\". This paramter is used in COMPUTE_LFC.')
        string(name: 'CONTROL_BARCODES', defaultValue: 'h-b', description: 'Type of control barcode ladder to be used in the pipeline. This defaults to \"h-b\".')

        // Parameters that we don't expect users to change
        separator(
          name: "Group_1",
          sectionHeader: "Do Not Edit",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleRed
        )

        // pipeline version
        string(name: 'GIT_BRANCH', defaultValue: 'main', description: 'Pipeline branch to use')
        booleanParam(name: 'USE_LATEST', defaultValue: true, description: 'Check this to use the most up to date version from the specified branch. If not checked, will use the specified commit.')
        string(name: 'COMMIT_ID', defaultValue: '', description: 'Specific commit ID to use (leave empty if using the latest commit in the branch or if already specified in the config file.)')

        // Sushi Input files
        string(name: 'RAW_COUNTS_UNCOLLAPSED', defaultValue: 'raw_counts_uncollapsed.csv', description: 'File name in BUILD_DIR containing the uncollapsed raw counts. This should be the file generated by Nori.')
        string(name: 'SAMPLE_META', defaultValue: 'sample_meta.csv', description: 'File name in BUILD_DIR of the sample meta.')
        string(name: 'CELL_SET_AND_POOL_META', defaultValue: 'cell_set_and_pool_meta.csv', description: 'Cell set and pool information for this run.')
        string(name: 'CELL_LINE_META', defaultValue: 'cell_line_meta.csv', description: 'File in BUILD_DIR containing cell line metadata')
        string(name: 'CONTROL_BARCODE_META', defaultValue: 'CB_meta.csv', description: 'Metadata for control barcodes')

        // Additional parameters ordered by when they first appear
        string(name: 'BARCODE_COL', defaultValue: 'forward_read_barcode', description: 'Name of the column containing the barcode sequence. The column containing the barcode sequence should have the same name across the Nori output file, the cell line metadata, and the CB metadata. This defaults to \"forward_read_barcode\", and the paramter is first used in COLLATE_FASTQ_READS.')
        string(name: 'SEQUENCING_INDEX_COLS', defaultValue: 'flowcell_names,index_1,index_2', description: 'List of sequencing related columns found in the sample meta. These columns are used to uniquely identify every PCR well in the run. Default value is \"flowcell_names,index_1,index_2\", and the parameter is used in COLLATE_FASTQ_READS.')
        string(name: 'ID_COLS', defaultValue: 'pcr_plate,pcr_well', description: 'List of columns found in the sample meta that are used to create a unique ID for each sample-replicate. This defaults to \"pcr_plate,pcr_well\", but could be any combination sample meta columns that uniquely identifies every well in the run. This parameter is first used in COLLATE_FASTQ_READS.')
        string(name: 'CHUNK_SIZE', defaultValue: '10000000', description: 'Number of rows for a chunk. Due to the large size of the Nori output, some actions are performed in chunks to conserve memory. This parameter sets the size of a chunk and defaults to 10^6 or \"10000000\". This paramter is first used in COLLATE_FASTQ_READS.')
        string(name: 'LOW_ABUNDANCE_THRESHOLD', defaultValue: '20', description: 'Threshold for unknown barcodes. Unknown barcodes below this threshold will be deidentified and their counts will be included under the term unknown_low_abundance_barcode. This paramter defaults to \"20\" and is used in COLLATE_FASTQ_READS.')
        string(name: 'PSEUDOCOUNT', defaultValue: '20', description: 'Pseudocount value added to all reads before log transformations. This defaults to \"20\" and is used in CBNORMALIZE.')
        string(name: 'CELL_LINE_COLS', defaultValue: 'pool_id,depmap_id,ccle_name', description: 'List of columns across the metadata files that are used to identify a unique cell line. This defaults to \"pool_id,depmap_id,ccle_name\", but can also include \"cell_set\" or descriptive columns like \"project_code\" that you would like to pass through the pipeline. This parameter is first used in COMPUTE_LFC.')
        string(name: 'COUNT_COL_NAME', defaultValue: 'normalized_n', description: 'Name of the numerical column that should be used to compute log2 fold change values. This defaults to \"normalized_n\" and is used in COMPUTE_LFC.')
        string(name: 'COUNT_THRESHOLD', defaultValue: '40', description: 'Threshold for filtering the negative controls. In the negative control conditions, cell lines whose median counts are below this threshold are not confidently detected and thus are dropped. This defaults to \"40\" and is used in COMPUTE_LFC.')

        // Files created by sushi
        string(name: 'PRISM_BARCODE_COUNTS', defaultValue: 'prism_barcode_counts.csv', description: 'Filename in BUILD_DIR containing PRISM barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'UNKNOWN_BARCODE_COUNTS', defaultValue: 'unknown_barcode_counts.csv', description: 'Filename in BUILD_DIR containing unknown barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'ANNOTATED_COUNTS', defaultValue: 'annotated_counts.csv', description: 'Filename in BUILD_DIR containing annotated counts. This file is creaed by FILTER_COUNTS.')
        string(name: 'FILTERED_COUNTS', defaultValue: 'filtered_counts.csv', description: 'Filename in BUILD_DIR containing filtered counts. This file is created by FILTER_COUNTS.')
        string(name: 'NORMALIZED_COUNTS', defaultValue: 'normalized_counts.csv', description: 'Filename in BUILD_DIR containing normalized counts. This file is created by CBNORMALIZE.')
        string(name: 'LFC', defaultValue: 'l2fc.csv', description: 'Filename containing log2 fold change values. This file is created by COMPUTE_LFC.')
        string(name: 'COLLAPSED_LFC', defaultValue: 'collapsed_l2fc.csv', description: 'Filename in BUILD_DIR containing replicate collapsed l2fc values. This file is created by COLLAPSED_LFC.')
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
                        'SEQ_TYPE', 'API_URL', 'BUILD_DIR', 'INDEX_1', 'INDEX_2', 'BARCODE_SUFFIX', 'CREATE_CELLDB_METADATA',
                        'BUILD_NAME', 'CONVERT_SUSHI', 'RUN_EPS_QC', 'REMOVE_DATA', 'DAYS',
                        'COUNTS', 'SCREEN', 'CONTROL_BARCODES',

                        // sushi input files
                        'RAW_COUNTS_UNCOLLAPSED', 'SAMPLE_META', 'CELL_SET_AND_POOL_META', 'CELL_LINE_META', 'CONTROL_BARCODE_META',

                        // sushi output files
                        'PRISM_BARCODE_COUNTS', 'UNKNOWN_BARCODE_COUNTS', 'ANNOTATED_COUNTS', 'FILTERED_COUNTS', 'NORMALIZED_COUNTS', 
                        'LFC', 'COLLAPSED_LFC',

                        // collate_fastq_reads parameters
                        'SEQUENCING_INDEX_COLS', 'ID_COLS', 'BARCODE_COL', 'LOW_ABUNDANCE_THRESHOLD', 'CHUNK_SIZE', 'REVERSE_INDEX2',

                        // normalize parameters
                        'PSEUDOCOUNT',

                        // compute_l2fc paramters
                        'SIG_COLS', 'CONTROL_COLS', 'CELL_LINE_COLS', 'COUNT_COL_NAME', 'CTL_TYPES', 'COUNT_THRESHOLD'
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
                        if (params.CREATE_SAMPLE_META) {
                            scriptsToRun.add('create_sample_meta.sh')
                        }
                        if (params.CREATE_CELLDB_METADATA) {
                            scriptsToRun.add('create_celldb_metadata.sh')
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
