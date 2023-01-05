
init:
	pip install -r requirements.txt

upgrade-dependencies:
	rm -f requirements.txt
	pip freeze > requirements.txt


