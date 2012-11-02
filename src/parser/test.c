#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bookshelf_IO.h"

int main(int argc, char *argv[]) {
		char auxFile[BUFFERSIZE], benchmarkPath[BUFFERSIZE], placefile[BUFFERSIZE];

  if(argc != 4) {
      printf("Usage: %s <benchmark_dir> <aux_file> <placement_file>\n",
             argv[0]);
      printf("    <benchmark_dir> is the benchmark file directory.\n");
      printf("    <aux_file> is the bookshelf format auxiliary file.\n");
      printf("    <placement_file> is the placement file.\n");
      exit(1);
  }    

  strcpy(benchmarkPath, argv[1]);
  strcpy(auxFile, argv[2]);
  strcpy(placefile, argv[3]);

  readAuxFile(benchmarkPath, auxFile);
  createHash(benchmarkPath, nodesFile);
  readNodesFile(benchmarkPath, nodesFile);
  readNetsFile(benchmarkPath, netsFile);
  readPlFile(benchmarkPath, placefile);
  freeHash();

	std::cout << movableNodes << " " << numTerminals << std::endl;
	std::cout << cellName[1] << "!\n";
	std::cout << cellName[movableNodes+numTerminals] << "!\n";
	for (int i = 0; i <= movableNodes+numTerminals; i++)
		std::cout << cellName[i] << "! " << i <<std::endl;
}
