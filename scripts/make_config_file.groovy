import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

String sectionHeaderStyleGreen = ' color: white; background: green; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 30px '
String sectionHeaderStyleBlue = ' color: white; background: blue; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 15px'
String sectionHeaderStyleRed = ' color: white; background: red; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 30px'
String separatorStyleCss = ' border: 0; border-bottom: 1px dashed #ccc; background: #999; '

pipeline {
    agent any

    parameters {
        separator(
            name: "analytics_modules",
            sectionHeader: "Analytics Modules",
            separatorStyle: separatorStyleCss,
            sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'DRC', defaultValue: false, description: 'Generate dose response curves.')
        booleanParam(name: 'UNIVARIATE_BIOMARKER', defaultValue: false, description: 'Run univariate biomarker analysis.')
        booleanParam(name: 'MULTIVARIATE_BIOMARKER', defaultValue: false, description: 'Run multivariate biomarker analysis.')
        booleanParam(name: 'LFC_BIOMARKER', defaultValue: false, description: 'Run LFC biomarker analysis.')
        booleanParam(name: 'AUC_BIOMARKER', defaultValue: false, description: 'Run AUC biomarker analysis.')
        booleanParam(name: 'QC_IMAGES', defaultValue: true, description: 'Check this to trigger the QC images job.')

        separator(
            name: "portal_prep",
            sectionHeader: "Portal Prep",
            separatorStyle: separatorStyleCss,
            sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'CONVERT_SUSHI', defaultValue: false, description: 'Convert output column headers to format for MTS pipeline and upload to s3.')
        string(name: 'DAYS', defaultValue: '', description: 'If running the sushi_to_mts module, provide any days/timepoints (separated by commas) that should be dropped from output data. No quotes needed (ie, 2,8).')

        separator(
            name: "metadata",
            sectionHeader: "Metadata",
            separatorStyle: separatorStyleCss,
            sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'CREATE_CELLDB_METADATA', defaultValue: true, description: 'Get cell line and pool metadata from cellDB.')
        booleanParam(name: 'CREATE_SAMPLE_META', defaultValue: false, description: 'Get sample metadata from COMET, project must be registered and all metadata steps completed.')

        // Active Choice Parameter for dynamically fetched screens
        activeChoiceParam(
            name: 'SCREEN',
            description: 'Select the screen name from COMET. Dynamically fetched from the API.',
            choiceType: 'SINGLE_SELECT',
            script: [
                classpath: [],
                script: '''
                    import groovy.json.JsonSlurper
                    import java.net.URL
                    import java.net.HttpURLConnection

                    try {
                        // API endpoint
                        String apiUrl = 'https://api.clue.io/api/prism_screens'

                        // Make the HTTP GET request
                        URL url = new URL(apiUrl)
                        HttpURLConnection connection = (HttpURLConnection) url.openConnection()
                        connection.setRequestMethod('GET')

                        if (connection.responseCode == 200) {
                            String response = url.text
                            def json = new JsonSlurper().parseText(response)
                            return json.collect { it.name } // Extract screen names
                        } else {
                            return ["Error: Unable to fetch screens (HTTP ${connection.responseCode})"]
                        }
                    } catch (Exception e) {
                        return ["Error fetching screens: ${e.message}"]
                    }
                '''
            ]
        )

        separator(
            name: "build_details",
            sectionHeader: "Build Details",
            separatorStyle: separatorStyleCss,
            sectionHeaderStyle: sectionHeaderStyleBlue
        )
        string(name: 'BUILD_DIR', defaultValue: '/cmap/obelix/pod/prismSeq/', description: 'Output path to deposit build. Format should be /directory/PROJECT_CODE/BUILD_NAME')
        string(name: 'BUILD_NAME', defaultValue: '', description: 'Build name; used to name output files.')
        string(name: 'SEQ_TYPE', defaultValue: 'DRAGEN', description: 'Choose DRAGEN, MiSeq, HiSeq, or NovaSeq. MiSeq and HiSeq/NovaSeq return files named differently. This setting sets the INDEX_1, INDEX_2, and BARCODE_SUFFIX parameters in fastq2readcount. Select DRAGEN if fastq files are from the DRAGEN pipeline from GP. Choosing NovaSeq reverses index 2.')

        // Additional parameters ordered by when they first appear
        string(name: 'BARCODE_COL', defaultValue: 'forward_read_barcode', description: 'Name of the column containing the barcode sequence. The column containing the barcode sequence should have the same name across the Nori output file, the cell line metadata, and the CB metadata. This defaults to \"forward_read_barcode\", and the parameter is first used in COLLATE_FASTQ_READS.')
        string(name: 'SEQUENCING_INDEX_COLS', defaultValue: 'flowcell_names,index_1,index_2', description: 'List of sequencing related columns found in the sample meta. These columns are used to uniquely identify every PCR well in the run. Default value is \"flowcell_names,index_1,index_2\", and the parameter is used in COLLATE_FASTQ_READS.')
        string(name: 'ID_COLS', defaultValue: 'pcr_plate,pcr_well', description: 'List of columns found in the sample meta that are used to create a unique ID for each sample-replicate. This defaults to \"pcr_plate,pcr_well\", but could be any combination sample meta columns that uniquely identifies every well in the run. This parameter is first used in COLLATE_FASTQ_READS.')
        string(name: 'CHUNK_SIZE', defaultValue: '10000000', description: 'Number of rows for a chunk. Due to the large size of the Nori output, some actions are performed in chunks to conserve memory. This parameter sets the size of a chunk and defaults to 10^6 or \"10000000\". This parameter is first used in COLLATE_FASTQ_READS.')
        string(name: 'LOW_ABUNDANCE_THRESHOLD', defaultValue: '20', description: 'Threshold for unknown barcodes. Unknown barcodes below this threshold will be deidentified and their counts will be included under the term unknown_low_abundance_barcode. This parameter defaults to \"20\" and is used in COLLATE_FASTQ_READS.')
        string(name: 'PSEUDOCOUNT', defaultValue: '20', description: 'Pseudocount value added to all reads before log transformations. This defaults to \"20\" and is used in CBNORMALIZE.')
        string(name: 'CELL_LINE_COLS', defaultValue: 'pool_id,depmap_id,lua', description: 'List of columns across the metadata files that are used to identify a unique cell line. This defaults to \"pool_id,depmap_id,lua\", but can also include \"cell_set\" or descriptive columns like \"project_code\" that you would like to pass through the pipeline. This parameter is first used in COMPUTE_LFC.')
        string(name: 'COUNT_COL_NAME', defaultValue: 'normalized_n', description: 'Name of the numerical column that should be used to compute log2 fold change values. This defaults to \"normalized_n\" and is used in COMPUTE_LFC.')
        string(name: 'COUNT_THRESHOLD', defaultValue: '40', description: 'Threshold for filtering the negative controls. In the negative control conditions, cell lines whose median counts are below this threshold are not confidently detected and thus are dropped. This defaults to \"40\" and is used in COMPUTE_LFC.')
        string(name: 'L2FC_COLUMN', defaultValue: 'l2fc', description: 'Name of the column containing the log2 fold change values used in DRC. This defaults to \"l2fc\".')
        string(name: 'VIABILITY_CAP', defaultValue: '1.5', description: 'Cap for viability values used when computing LFC. This defaults to \"1.5\".')

        // Files created by sushi
        string(name: 'PRISM_BARCODE_COUNTS', defaultValue: 'prism_barcode_counts.csv', description: 'Filename in BUILD_DIR containing PRISM barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'UNKNOWN_BARCODE_COUNTS', defaultValue: 'unknown_barcode_counts.csv', description: 'Filename in BUILD_DIR containing unknown barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'ANNOTATED_COUNTS', defaultValue: 'annotated_counts.csv', description: 'Filename in BUILD_DIR containing annotated counts. This file is created by FILTER_COUNTS.')
        string(name: 'FILTERED_COUNTS', defaultValue: 'filtered_counts.csv', description: 'Filename in BUILD_DIR containing filtered counts. This file is created by FILTER_COUNTS.')
        string(name: 'NORMALIZED_COUNTS', defaultValue: 'normalized_counts.csv', description: 'Filename in BUILD_DIR containing normalized counts. This file is created by CBNORMALIZE.')
        string(name: 'LFC', defaultValue: 'l2fc.csv', description: 'Filename containing log2 fold change values. This file is created by COMPUTE_LFC.')
        string(name: 'COLLAPSED_LFC', defaultValue: 'collapsed_l2fc.csv', description: 'Filename in BUILD_DIR containing replicate collapsed l2fc values. This file is created by COLLAPSED_LFC.')
        string(name: 'API_URL', defaultValue: 'https://api.clue.io/api/', description: 'API URL')
    }

    environment {
        CONFIG_FILE_PATH = "${params.BUILD_DIR}/config.json"
    }

    stages {
        stage('Checkout') {
            steps {
                script {
                    echo "Checking out the repository..."
                }
            }
        }

        stage('Generate JSON Config') {
            steps {
                script {
                    def config = [:]

                    // Add parameters to the config
                    config['SCREEN'] = params.SCREEN
                    config['CREATE_SAMPLE_META'] = params.CREATE_SAMPLE_META

                    writeFile file: env.CONFIG_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                    echo "Configuration generated: ${config}"
                }
            }
        }

        stage('Run Scripts in Container') {
            steps {
                script {
                    if (params.CREATE_SAMPLE_META) {
                        echo "Running create_sample_meta script with screen: ${params.SCREEN}"
                        sh """
                            ./scripts/create_sample_meta.sh --screen ${params.SCREEN}
                        """
                    }
                }
            }
        }
    }

    post {
        success {
            echo "Pipeline completed successfully."
        }
        failure {
            echo "Pipeline failed."
        }
    }
}
