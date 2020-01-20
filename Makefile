THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))
IMAGE = euclid:1.0
CONTAINER = euclid-container

build: 
	docker build --tag $(IMAGE)  \
	    -f Dockerfile \
        .

run:
	docker run --rm \
		-it \
		-p 7774:22 \
	    -v $(THIS_FILE_DIR):/opt/project \
	    --name $(CONTAINER) \
	    $(IMAGE)
