import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

pipeline {
    agent any

    // Define parameters that can be edited via the Jenkins UI
    parameters {
        string(name: 'BUILD_DIR', defaultValue: '', description: 'Output path to deposit build. Format should be /directory/PROJECT_CODE/BUILD_NAME')
        booleanParam(name: 'TRIGGER_BUILD', defaultValue: true, description: 'Check this to trigger the build.')
    }

    environment {
        CONFIG_FILE_PATH = "${params.BUILD_DIR}/config.json"
    }

    stages {
        stage('Check Config File') {
            steps {
                script {
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        def config = new JsonSlurper().parseText(configText)
                        config.each { key, value ->
                            env.setProperty(key, value.toString())
                        }
                        echo "Loaded config from ${env.CONFIG_FILE_PATH}"
                    } else {
                        def config = [
                            GIT_BRANCH: params.GIT_BRANCH,
                            COMMIT_HASH: params.COMMIT_HASH,
                            BUILD_DIR: params.BUILD_DIR,
                            BUILD_NAME: params.BUILD_NAME,
                            SEQ_TYPE: params.SEQ_TYPE,
                            CTL_TYPES: params.CTL_TYPES,
                            DAYS: params.DAYS,
                            CELL_SET_META: params.CELL_SET_META,
                            ID_COLS: params.ID_COLS,
                            SAMPLE_COLS: params.SAMPLE_COLS,
                            SIG_COLS: params.SIG_COLS,
                            SEQUENCING_INDEX_COLS: params.SEQUENCING_INDEX_COLS,
                            CONTROL_BARCODE_META: params.CONTROL_BARCODE_META,
                            CONTROL_COLS: params.CONTROL_COLS,
                            REMOVE_DATA: params.REMOVE_DATA,
                            COUNT_COL_NAME: params.COUNT_COL_NAME,
                            RUN_NORM: params.RUN_NORM,
                            PULL_POOL_ID: params.PULL_POOL_ID,
                            CONVERT_SUSHI: params.CONVERT_SUSHI,
                            RUN_EPS_QC: params.RUN_EPS_QC,
                            SAMPLE_META: params.SAMPLE_META,
                            COUNT_THRESHOLD: params.COUNT_THRESHOLD,
                            PSEUDOCOUNT: params.PSEUDOCOUNT,
                            CELL_LINE_META: params.CELL_LINE_META,
                            RAW_COUNTS: params.RAW_COUNTS,
                            FILTERED_COUNTS: params.FILTERED_COUNTS,
                            LFC: params.LFC,
                            ANNOTATED_COUNTS: params.ANNOTATED_COUNTS,
                            NORMALIZED_COUNTS: params.NORMALIZED_COUNTS,
                            COLLAPSED_VALUES: params.COLLAPSED_VALUES,
                            API_URL: params.API_URL,
                            API_KEY: sh(script: 'cat /local/jenkins/.clue_api_key', returnStdout: true).trim()
                        ]
                        writeFile file: env.CONFIG_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                        echo "Generated config.json: ${config}"
                    }
                }
            }
        }

        stage('Checkout') {
            steps {
                script {
                    def branchSpec = env.COMMIT_HASH?.trim() ? env.COMMIT_HASH : "*/${env.GIT_BRANCH}"
                    checkout([
                        $class: 'GitSCM',
                        branches: [[name: branchSpec]],
                        userRemoteConfigs: scm.userRemoteConfigs
                    ])
                }
            }
        }

        stage('Add Commit Hash to Config') {
            steps {
                script {
                    def commitHash = sh(script: 'git rev-parse HEAD', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg commit "$commitHash" '. + {COMMIT: \$commit}' $configFilePath > ${configFilePath}.tmp && mv ${configFilePath}.tmp $configFilePath
                    """
                    echo "Added commit hash to config.json: $commitHash"
                }
            }
        }
    }

    post {
        success {
            script {
                if (env.TRIGGER_BUILD.toBoolean()) {
                    def configText = readFile(file: env.CONFIG_FILE_PATH)
                    def config = new JsonSlurper().parseText(configText)
                    def nextJobParams = config.collect { key, value ->
                        if (value == 'true' || value == 'false') {
                            booleanParam(name: key, value: value.toBoolean())
                        } else {
                            string(name: key, value: value.toString())
                        }
                    }
                    build job: 'create_celldb_metadata_podman', wait: false, parameters: nextJobParams
                } else {
                    echo 'Next build not triggered because TRIGGER_BUILD is false.'
                }
            }
        }
    }
}
