nextflow.enable.dsl=2

def getServerPopulations(tbm_gen_server) {
    def serverPopulations = [
        "s115": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s107": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature", "test_chr12"],
        "s66": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s203": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s204": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s243": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature", "w115", "w66", "w203", "w204", "w243"]
    ]

    return serverPopulations[tbm_gen_server]
}

def getTigerJarConfig(tiger_jar_name, home_dir, app_name_override = null) {
    def tiger_jar_versions = [
        "TIGER_F3_2M_scan_20260129.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20251121.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20251118.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20251016.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20250915.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_20250526.jar": [
            java_version: "java8", 
            app_name: "FastCall2"
        ],
        "TIGER_PD_20260112.jar": [
            java_version: "java17", 
            app_name: "PopDep"
        ],
        "TIGER_PD_20260130.jar": [
            java_version: "java17", 
            app_name: "PopDep"
        ],
        "TIGER_PD_20260615.jar": [
            java_version: "java17",
            app_name: "PopDepFull"
        ]
    ]

    def config = tiger_jar_versions[tiger_jar_name]
    if (!config) {
        throw new IllegalArgumentException("Unknown TIGER jar: ${tiger_jar_name}")
    }
    def tigerJar = file("${home_dir}/lib/${tiger_jar_name}")
    
    return [
        path: tigerJar,
        app_name: app_name_override ?: config.app_name,
        java_version: config.java_version
    ]
}
