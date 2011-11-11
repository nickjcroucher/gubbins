newicktree: seqMain.o seqUtil.o Newickform.o
	gcc -o newicktree seqMain.o seqUtil.o Newickform.o
seqMain.o: seqMain.c
	gcc -c seqMain.c
seqUtil.o: seqUtil.c
	gcc -c seqUtil.c
Newickform.o: Newickform.c
	gcc -c Newickform.c

clean:
	-rm *.o newicktree
