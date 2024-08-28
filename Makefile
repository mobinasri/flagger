# Definitions
repository = mobinasri
identifier = flagger
version = v1.0_alpha
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	docker tag ${name}:${tag} ${name}:${version}

push: build
	docker push ${name}:${tag}
	docker push ${name}:latest
	docker push ${name}:${version}
