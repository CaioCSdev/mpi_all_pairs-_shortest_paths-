.PHONY: run run-cluster deploy

fox: main.c main.h
	mpicc -Wall --pedantic -g -o fox main.c

run: fox
	@mpirun fox

deploy: fox
	rsync --rsh 'ssh -J up202410254@ssh.alunos.dcc.fc.up.pt' --human-readable --recursive . up202410254@L102:/net/areas/homes/up202410254/classes/parallel/pa1

run-cluster: deploy
	@ssh -J up202410254@ssh.alunos.dcc.fc.up.pt up202410254@L102 -t 'cd classes/parallel/pa1; make fox ; mpirun --hostfile hostnames --use-hwthread-cpus fox'

