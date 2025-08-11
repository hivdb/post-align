requirements.txt: Pipfile.lock
	@pipenv requirements --from-pipfile > requirements.txt
@sed -i.bak -- '/^-e \.$/d' requirements.txt && rm -f requirements.txt.bak

test-unit:
	@pytest tests/unit

test-component:
	@behave tests/component

test: test-unit test-component

build-docker-builder: requirements.txt
	@docker pull ubuntu:18.04
	# @docker build . --no-cache -t hivdb/post-align-builder:latest
	@docker build . -t hivdb/post-align-builder:latest

push-docker-builder: build-docker-builder
	@docker push hivdb/post-align-builder:latest

inspect-builder:
	@docker run \
		--mount type=bind,source=$(PWD),target=/app/post-align \
		--workdir /app --rm -it \
		hivdb/post-align-builder:latest \
		bash

dist/linux-amd64: $(shell find postalign -type f -path "*.py" | sed 's#\([| ]\)#\\\1#g')
	@docker run \
		--mount type=bind,source=$(PWD),target=/app/post-align \
		--workdir /app --rm -it \
		hivdb/post-align-builder:latest \
		/build.sh

dist/postalign_linux-amd64.tar.gz: dist/linux-amd64 
	@cd dist/linux-amd64 && \
		rm -f mpostalign_linux-amd64.tar.gz && \
		tar zcf mpostalign_linux-amd64.tar.gz mpostalign && \
		mv mpostalign_linux-amd64.tar.gz ..

dist: dist/postalign_linux-amd64.tar.gz
build: dist

.PHONY: test-unit test-component test build-docker-builder push-docker-builder
