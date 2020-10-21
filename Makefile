requirement.txt: Pipfile.lock
	@pipenv lock --dev --requirements > requirements.txt

build-docker-builder: requirement.txt
	@docker pull ubuntu:18.04
	@docker build . --no-cache -t hivdb/post-align-builder:latest

push-docker-builder: build-docker-builder
	@docker push hivdb/post-align-builder:latest

dist/linux-amd64: $(shell find postalign -type f -path "*.py" | sed 's#\([| ]\)#\\\1#g')
	@rm -rf dist/linux-amd64
	@mkdir -p dist/linux-amd64
	@docker run \
		--mount type=bind,source=$(PWD)/postalign,target=/app/postalign \
		--mount type=bind,source=$(PWD)/dist/linux-amd64,target=/app/dist \
		--workdir /app --rm -it \
		hivdb/post-align-builder:latest \
		pyinstaller postalign/entry.py -n postalign
	@cd dist/linux-amd64 && \
		tar zcf postalign_linux-amd64.tar.gz postalign && \
		mv postalign_linux-amd64.tar.gz ..

dist/darwin-amd64: $(shell find postalign -type f -path "*.py" | sed 's#\([| ]\)#\\\1#g')
	@test "$(shell uname)" = "Darwin" && \
		pipenv run pyinstaller postalign/entry.py -n postalign --distpath ./dist/darwin-amd64
	@rm ./postalign.spec
	@cd dist/darwin-amd64 && \
		tar zcf postalign_darwin-amd64.tar.gz postalign && \
		mv postalign_darwin-amd64.tar.gz ..

dist/postalign_linux-amd64.tar.gz: dist/linux-amd64 dist/darwin-amd64

dist: dist/postalign_linux-amd64.tar.gz
build: dist

.PHONY: build-builder