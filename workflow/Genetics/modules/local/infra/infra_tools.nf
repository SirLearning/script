nextflow.enable.dsl=2

include { getTigerJarConfig } from './infra_tiger.nf'

def getJavaSetupScript(javaVersion, javaLibDir) {
    def javaVersionMap = [
        "java8": "jdk-8",
        "java11": "jdk-11", 
        "java17": "jdk-17",
        "java21": "jdk-21"
    ]
    
    def javaDir = javaVersionMap[javaVersion]
    if (!javaDir) {
        def validVersions = javaVersionMap.keySet().join(", ")
        throw new Exception("Unknown Java version: ${javaVersion}. Valid options: ${validVersions}")
    }
    
	def libDir = javaLibDir ?: "${params.user_dir}/lib/jvm"
    def javaHome = "${libDir}/${javaDir}"
    
    return """
    # Switch to ${javaVersion} (${javaDir})
    export JAVA_HOME=${javaHome}
    export PATH=\$JAVA_HOME/bin:\$PATH
    
    # Verify Java version
    if [ ! -d "\$JAVA_HOME" ]; then
        echo "Error: Java directory not found: \$JAVA_HOME"
        echo "Available Java versions in ${javaLibDir}:"
        ls -la ${javaLibDir}/ || echo "Java lib directory not accessible"
        exit 1
    fi
    
    java -version 2>&1 | head -3
    """.stripIndent()
}
def getSoftwareConfig(home_dir, user_dir, tiger_jar, workshop_jar, samtools) {

    def java_lib_dir = "${user_dir}/lib/jvm"

    def resolved_workshop_jar = "${home_dir}/lib/${workshop_jar}"

    def tiger_jar_config = getTigerJarConfig(tiger_jar, home_dir)
    def workshop_jar_config = getWorkshopJarConfig(resolved_workshop_jar, java_lib_dir)
    def samtools_config = getSamtoolsConfig(samtools)
    
    return [
        tiger_jar_config: tiger_jar_config,
        workshop_jar_config: workshop_jar_config,
        samtools_config: samtools_config
    ]
}

def getWorkshopJarConfig(workshopJarPath, java_lib_dir) {
    return [
        path: workshopJarPath,
        java_version: "java21",
        java_lib_dir: java_lib_dir
    ]
}

def getSamtoolsConfig(samtoolsPath) {
    return [
        path: samtoolsPath
    ]
}
def validatePaths(pathMap) {
    def errors = []
    def isValid = true
    
    pathMap.each { name, path ->
        if (!path) {
            errors << "${name} is not specified"
            isValid = false
        } else {
            def pathFile = file(path)
            if (!pathFile.exists()) {
                errors << "${name} does not exist: ${path}"
                isValid = false
            }
        }
    }
    
    return [isValid: isValid, errors: errors]
}
