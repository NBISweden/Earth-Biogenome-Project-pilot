#!/bin/bash

IMAGE_VERSION=1.3

IMAGE_NAME="ghcr.io/nbisweden/fastk_genescopefk_merquryfk"

# Build

	docker build --network=host --platform=linux/amd64 -f ./fastk_genescopefk_merquryfk.Dockerfile -t $IMAGE_NAME:$IMAGE_VERSION . || {
			echo "Error building image, exited script"
			exit 1
	}

# Push

	# docker push $IMAGE_NAME:$IMAGE_VERSION || {
	#         echo "Error pushing image, exited script"
	#         exit 1
	# }

# Debug build

	# docker build --network=host --platform=linux/amd64 -f ./fastk_genescopefk_merquryfk.Dockerfile -t $IMAGE_NAME:$IMAGE_VERSION . > build.log 2>&1 || {
	#       echo "Error building image, exited script"
	#       exit 1
	# }
