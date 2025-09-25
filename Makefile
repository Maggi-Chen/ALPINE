NAME=ALPINE
VERSION=$(shell git branch --show-current)
git_hash=$(shell git rev-parse --short HEAD)

LAUNCHER_NAME=transgene-launcher

REGISTRY= # Add AWS registry here

# ARVADOS VARIABLES
MASTER_PROJECT_UUID= # Add UUID
LAUNCHER_PIPELINE_UUID= # Add UUID
PER_SAMPLE_CLASSIFICATION_WORKFLOW= # Add UUID
MERGE_CLASSIFICATION_WORKFLOW= # Add UUID


# SEVENBRIDGES VARIABLES
SBG_MASTER_PROJECT= # Add SBG ID
SBG_REGISTRY= # Add SBG registry 

all: help 
.PHONY: help login build push test

help:
	@echo "Usage: make [command]"
	@echo "  Command      Description"
	@echo "  -------      -----------"
	@echo "  login        Log in to ECR"
	@echo "  build        Build Docker image"
	@echo "  push         Push Docker image to ECR"
	@echo " "
	@echo " sevenbridges-build-launcher  Build launcher Docker image and CWL for SevenBridges"
	@echo " sevenbridges-push-launcher  Push launcher Docker image and CWL for SevenBridges"
	@echo " "
	@echo " arvados-build-launcher  Build launcher Docker image for Arvados"
	@echo " arvados-push-launcher   Push launcher Docker image for Arvados"
	@echo " arvados-push-workflows  Push workflow for WES pipeline to Arvados"

login:
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${REGISTRY}

build: 
	cd Dockerfiles; ./build_images.sh; cd ..

test:
	cd launcher/arvados && pytest --cov-config=.coveragerc --cov src/ && coverage html -i

push: login
	##cd Dockerfiles; ./deploy_images.sh; cd ..

build-launcher:
	PYTHONPATH="launcher/src/" pylint --max-line-length=120 --recursive=true launcher/src/ --exit-zero
	docker build --platform linux/amd64 -t $(LAUNCHER_NAME):$(VERSION) launcher/
	sed 's/__NAME__/$(LAUNCHER_NAME)/g' launcher/launcher.cwl | sed 's/__VERSION__/$(VERSION)/' > launcher/launcher.$(VERSION).cwl

sevenbridges-push-launcher: build-launcher
	docker login images.sbgenomics.com
	docker tag $(LAUNCHER_NAME):$(VERSION) $(SBG_REGISTRY)/$(LAUNCHER_NAME):$(VERSION)
	docker push $(SBG_REGISTRY)/$(LAUNCHER_NAME):$(VERSION)
	sbpack default $(SBG_MASTER_PROJECT)/$(LAUNCHER_NAME) launcher/launcher.$(VERSION).cwl

sevenbridges-push-workflow:
	sbpack default $(SBG_MASTER_PROJECT)/per-sample-classification CWL/aav_per_sample_workflow.cwl
	sbpack default $(SBG_MASTER_PROJECT)/merge-classification-table CWL/merge_count_table.cwl

arvados-push-launcher: build-launcher
	docker login images.sbgenomics.com
	docker tag $(LAUNCHER_NAME):$(VERSION) $(SBG_REGISTRY)/$(LAUNCHER_NAME):$(VERSION)
	docker push $(SBG_REGISTRY)/$(LAUNCHER_NAME):$(VERSION)
	arvados-cwl-runner --update-workflow $(LAUNCHER_PIPELINE_UUID) --project-uuid $(MASTER_PROJECT_UUID) launcher/launcher.$(VERSION).cwl

arvados-push-workflows:
	arvados-cwl-runner --update-workflow $(PER_SAMPLE_CLASSIFICATION_WORKFLOW) CWL/aav_per_sample_workflow.cwl
	arvados-cwl-runner --update-workflow $(MERGE_CLASSIFICATION_WORKFLOW) CWL/merge_count_table.cwl


