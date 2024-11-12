.PHONY: run run-cluster deploy

fox: main.c main.h
	mpicc -Wall --pedantic -g -o fox main.c -lm

run: fox
	@mpirun --use-hwthread-cpus fox input

deploy: fox
	rsync --rsh 'ssh -J up202410254@ssh.alunos.dcc.fc.up.pt' --human-readable --recursive . up202410254@L102:/net/areas/homes/up202410254/classes/parallel/pa1

run-cluster: deploy
	@ssh -J up202410254@ssh.alunos.dcc.fc.up.pt up202410254@L102 -t 'cd classes/parallel/pa1; make fox ; mpirun --hostfile hostnames --use-hwthread-cpus fox input'

test: fox
	@mpirun fox input > code.output
	@diff code.output output -q
