all:
	@echo '"tosvg" install start!'
	npm install
	npm link
	@echo ''
	@echo '*************************************************'
	@echo 'congratulations!'
	@echo '*************************************************'
	@echo '"tosvg" has been installed'
	@echo 'If you are using "nodenv" you may need to type "nodenv rehash" command'
	@echo ''
	@echo '[Usage]'
	@echo 'tosvg img.png img.svg'