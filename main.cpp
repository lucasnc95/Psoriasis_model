#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "OpenCLWrapper.h"
#include <iostream>
#include <stdio.h>


#define CELULAS_SINGLE_CELL_SIZEOF 8


#define CELULAS_SIZEOF 40
#define MALHA_TOTAL_CELULAS 5
// Offsets para as células na malha
#define CELULAS_L_OFFSET 0    // Células T - L
#define CELULAS_M_OFFSET 8    // Células dendríticas - M
#define CELULAS_S_OFFSET 16   // Citocinas - S
#define CELULAS_W_OFFSET 24   // Medicamento imunológico - W
#define CELULAS_K_OFFSET 32   // Queratinócitos - K


#define CELULA_L  0   // Células T (L)
#define CELULA_M  1   // Células dendríticas (M)
#define CELULA_S  2   // Citocinas (S)
#define CELULA_W  3   // Medicamento imunológico (W)
#define CELULA_K  4   // Queratinócitos (K)

//Informacoes de acesso à estrutura "parametrosMalha".
#define OFFSET_COMPUTACAO               0
#define LENGTH_COMPUTACAO               1
#define COMPRIMENTO_GLOBAL_X            2
#define COMPRIMENTO_GLOBAL_Y            3
#define COMPRIMENTO_GLOBAL_Z            4
#define MALHA_DIMENSAO_POSICAO_Z        5
#define MALHA_DIMENSAO_POSICAO_Y        6
#define MALHA_DIMENSAO_POSICAO_X        7
#define MALHA_DIMENSAO_CELULAS          8
#define NUMERO_PARAMETROS_MALHA         9



// Definição de semente fixa para reprodutibilidade
#define SEED 12345





void InicializarParametrosMalhaHIS(int *parametrosMalha,  int offsetComputacao,  int lengthComputacao,  int xMalhaLength,  int yMalhaLength,  int zMalhaLength)
{
	//parametrosMalha = new int[NUMERO_PARAMETROS_MALHA];

	(parametrosMalha)[OFFSET_COMPUTACAO] = offsetComputacao;
	(parametrosMalha)[LENGTH_COMPUTACAO] = lengthComputacao;
	(parametrosMalha)[COMPRIMENTO_GLOBAL_X] = xMalhaLength;
	(parametrosMalha)[COMPRIMENTO_GLOBAL_Y] = yMalhaLength;
	(parametrosMalha)[COMPRIMENTO_GLOBAL_Z] = zMalhaLength;
	(parametrosMalha)[MALHA_DIMENSAO_POSICAO_Z] = yMalhaLength*xMalhaLength*MALHA_TOTAL_CELULAS;
	(parametrosMalha)[MALHA_DIMENSAO_POSICAO_Y] = xMalhaLength*MALHA_TOTAL_CELULAS;
	(parametrosMalha)[MALHA_DIMENSAO_POSICAO_X] = MALHA_TOTAL_CELULAS;
	(parametrosMalha)[MALHA_DIMENSAO_CELULAS] = 1;
}




// Função para inicializar os pontos da malha para psoríase
void InicializarPontosPsoriase(float *malha, int *parametrosMalha, int K, int W, int L, int M, int S) {
    // Limpa a malha
    for(unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++) {
        for(unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++) {
            for(unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++) {
                (malha)[(CELULA_L * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
                (malha)[(CELULA_M * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
                (malha)[(CELULA_S * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
                (malha)[(CELULA_W * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
                (malha)[(CELULA_K * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
            }
        }
    }

    // Inicializa K (Queratinócitos)
    int centerX = parametrosMalha[COMPRIMENTO_GLOBAL_X] / 2;
    int centerY = parametrosMalha[COMPRIMENTO_GLOBAL_Y] / 2;
    int centerZ = parametrosMalha[COMPRIMENTO_GLOBAL_Z] / 2;

    // Preencher K no centro da malha
    if (K > 0) {
        for (int i = 0; i < K; i++) {
            int offsetX = (i % 3) - 1; // Adjacências em X
            int offsetY = (i / 3) % 3 - 1; // Adjacências em Y
            int offsetZ = i / 9 - 1; // Adjacências em Z

            int posX = centerX + offsetX;
            int posY = centerY + offsetY;
            int posZ = centerZ + offsetZ;

            if (posX >= 0 && posX < parametrosMalha[COMPRIMENTO_GLOBAL_X] &&
                posY >= 0 && posY < parametrosMalha[COMPRIMENTO_GLOBAL_Y] &&
                posZ >= 0 && posZ < parametrosMalha[COMPRIMENTO_GLOBAL_Z]) {
                (malha)[(CELULA_K * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (posZ * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (posY * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (posX * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
            }
        }
    }

    // Inicializa W (Medicamento imunológico) nas adjacências do centro
    if (W > 0) {
        for (int i = 0; i < W; i++) {
            int offsetX = (i % 3) - 1; // Adjacências em X
            int offsetY = (i / 3) % 3 - 1; // Adjacências em Y
            int offsetZ = i / 9 - 1; // Adjacências em Z

            int posX = centerX + offsetX;
            int posY = centerY + offsetY;
            int posZ = centerZ + offsetZ;

            if (posX >= 0 && posX < parametrosMalha[COMPRIMENTO_GLOBAL_X] &&
                posY >= 0 && posY < parametrosMalha[COMPRIMENTO_GLOBAL_Y] &&
                posZ >= 0 && posZ < parametrosMalha[COMPRIMENTO_GLOBAL_Z]) {
                (malha)[(CELULA_W * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (posZ * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (posY * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (posX * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
            }
        }
    }

    // Inicializa L, M e S em posições aleatórias
    for (int i = 0; i < L; i++) {
        int posX = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_X];
        int posY = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Y];
        int posZ = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Z];

        (malha)[(CELULA_L * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (posZ * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (posY * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (posX * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
    }

    for (int i = 0; i < M; i++) {
        int posX = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_X];
        int posY = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Y];
        int posZ = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Z];

        (malha)[(CELULA_M * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (posZ * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (posY * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (posX * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
    }

    for (int i = 0; i < S; i++) {
        int posX = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_X];
        int posY = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Y];
        int posZ = rand() % parametrosMalha[COMPRIMENTO_GLOBAL_Z];

        (malha)[(CELULA_S * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (posZ * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (posY * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (posX * parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
    }
}



void SaveFigure(float *malha, int *parametrosMalha, int xMalhaLength, int yMalhaLength, int zMalhaLength, int time)
{
    FILE *file;
    char filename[50];
    char charTime[10];
    
    // Criação do nome do arquivo com o tempo
    snprintf(filename, sizeof(filename), "results%d.vtk", time);
    
    // Abre o arquivo para escrita
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Não foi possível abrir o arquivo");
        return;
    }

    // Cabeçalho do arquivo VTK
    fprintf(file, "# vtk DataFile Version 2.0\n");
    fprintf(file, "Really cool data\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_GRID\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", xMalhaLength, yMalhaLength, zMalhaLength);
    fprintf(file, "POINTS %d float\n", xMalhaLength * yMalhaLength * zMalhaLength);

    // Escreve as coordenadas dos pontos
    for (unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++) {
        for (unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++) {
            for (unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++) {
                fprintf(file, "%d %d %d\n", x, y, z);
            }
        }
    }

    fprintf(file, "POINT_DATA %d\n", xMalhaLength * yMalhaLength * zMalhaLength);
    fprintf(file, "SCALARS volume_scalars float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    // Preenche os valores escalares
    
    for (unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++) {
        for (unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++) {
            for (unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++) {
                int index = (CELULA_L * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X]);
                
                if (index >= parametrosMalha[OFFSET_COMPUTACAO] * MALHA_TOTAL_CELULAS &&
                    index < (parametrosMalha[OFFSET_COMPUTACAO] + parametrosMalha[LENGTH_COMPUTACAO]) * MALHA_TOTAL_CELULAS) {
                    fprintf(file, "%f ", malha[index]);
                } else {
                    fprintf(file, "0.0 ");
                }
            }
            fprintf(file, "\n");
        }
    }

    fclose(file);
}


void LerPontosHIS(float *malha, int *parametrosMalha, int celulaDesejada)
{
    int counter = 0; // Contador para células impressas
    int gtz = 0;
    for (unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++)
    {
        for (unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++)
        {
            for (unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++)
            {
                // Calcula o índice da célula desejada
                int indice = (celulaDesejada * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + 
                             (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + 
                             (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + 
                             (x * parametrosMalha[MALHA_DIMENSAO_POSICAO_X]);

                // Verifica se o índice está no intervalo de computação
                if (indice >= parametrosMalha[OFFSET_COMPUTACAO] * MALHA_TOTAL_CELULAS && 
                    indice < (parametrosMalha[OFFSET_COMPUTACAO] + parametrosMalha[LENGTH_COMPUTACAO]) * MALHA_TOTAL_CELULAS)
                {
                    counter++; // Incrementa o contador
                    printf("%f ", malha[indice]); // Imprime o valor da célula desejada
                    if(malha[indice]>0)
                        gtz++;
                }
                else
                {
                    printf("%f ", 0.0f); // Imprime 0.0 se fora do intervalo
                }
            }
            printf("\n"); // Nova linha após cada camada z
        }
    }

    std::cout << "Células impressas: " << counter << std::endl; // Imprime o número total de células impressas
    std::cout << "Células impressas maiores que zero: " << gtz << std::endl; // Imprime o número total de células impressas
}


int main(int argc, char** argv) {
     
    OpenCLWrapper openCL(argc, argv);
    openCL.InitDevices("GPU_DEVICES", 10);  
    openCL.setKernel("kernel_psoriase.cl", "ProcessarPontos");
    
    int x = 100, y = 100, z = 100;
    int tam = x * y * z * MALHA_TOTAL_CELULAS;  // Tamanho correto da malha
    int *parametros = new int[NUMERO_PARAMETROS_MALHA];
    float *malha = new float[tam];  // Alocar a malha corretamente
    
    double	tempoInicio = MPI_Wtime();
    InicializarParametrosMalhaHIS(parametros, 0, (x * y * z), x, y, z);
    
   
    InicializarPontosPsoriase(malha, parametros, 100, 10, 30, 20, 25);
	int total_elements = x*y*z;
    
    
    // Configuração do balanceador de carga no OpenCLWrapper
    size_t sub = x * y  ;
    openCL.setLoadBalancer(sizeof(float), total_elements, MALHA_TOTAL_CELULAS, sub);
    
    // Alocar objetos de memória OpenCL
    int aMemObj = openCL.AllocateMemoryObject(NUMERO_PARAMETROS_MALHA * sizeof(int), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, parametros);
    int bMemObj = openCL.AllocateMemoryObject(tam * sizeof(float), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, malha);
    int cMemObj = openCL.AllocateMemoryObject(tam * sizeof(float), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, malha);
    
    float *malhaAux = new float[tam];  // Para armazenar os resultados
	
	
    // Definir atributos de kernel
	openCL.setAttribute(0, bMemObj);
    openCL.setAttribute(1, cMemObj);
	openCL.setAttribute(2, aMemObj);
	openCL.setBalancingTargetID(bMemObj);
	
    int contador = 0;
    for (int x = 0; x < 100000; x++) {
		
	
		 if (x % 2 == 0) {
            openCL.setAttribute(0, bMemObj);
            openCL.setAttribute(1, cMemObj);
            openCL.setBalancingTargetID(bMemObj);
        } 
		else {
            openCL.setAttribute(0, cMemObj);
            openCL.setAttribute(1, bMemObj);
            openCL.setBalancingTargetID(cMemObj);
        }

		if(x%1000  == 0){
        openCL.GatherResults(bMemObj, malhaAux);
        SaveFigure(malhaAux,parametros,10,10,10,contador);
        contador++;
		}
        
		openCL.ExecuteKernel();
		
    }

	openCL.GatherResults(bMemObj, malhaAux);
    LerPontosHIS(malhaAux, parametros, CELULA_L);
	  
	double tempoFim = MPI_Wtime();
	std::cout<<"Tempo execução:"<<tempoFim-tempoInicio<<std::endl;
    // Liberar memória alocada
    delete[] parametros;
    delete[] malha;
    delete[] malhaAux;
	
    return 0;
}
