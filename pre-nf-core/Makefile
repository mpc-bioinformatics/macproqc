docker-imgs:
	docker build -t ghcr.io/mpc-bioinformatics/mcquac:dev -f docker/Dockerfile .
	docker pull quay.io/medbioinf/comet-ms:v2024.01.0
	docker pull quay.io/biocontainers/pia:1.5.5--hdfd78af_0
	docker pull quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0
	docker pull mfreitas/tdf2mzml

comet-params:
	docker run --rm -it quay.io/medbioinf/comet-ms:v2024.01.0 bash -c 'comet -p > /dev/null; cat comet.params.new' > comet.params.new