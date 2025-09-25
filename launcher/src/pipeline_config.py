'''
Config settings specific to this pipeline
'''
PIPELINE = "Transgene"

__appname__ = f"{PIPELINE} Launcher"
__version__ = "0.1"

REFERENCE_PROJECT_NAME = f"{PIPELINE} Reference"

WORKFLOW_TO_OUTPUTS = {
    "per_sample": [
        ("count_table","Per_Sample"),
    ],
}

PIPELINE_PLATFORM_CONFIG = {
    "Arvados": {
        "per_sample_workflow": "", # Add Arvados workflow ID here
        "merge_workflow": "", # Add Arvados workflow ID here
    },
    "SevenBridges": {
        "per_sample_workflow": "", # Add SBG workflow ID here
        "merge_workflow": "", # Add SBG workflow ID here
        "api_endpoint": "https://", # Add SBG API endpoint here
        "token": "dummy",
    }
}
