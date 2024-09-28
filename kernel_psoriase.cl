// Tipos de células
#define MALHA_TOTAL_CELULAS   5
#define CELULA_L  0   // Células T (L)
#define CELULA_M  1   // Células dendríticas (M)
#define CELULA_S  2   // Citocinas (S)
#define CELULA_W  3   // Medicamento imunológico (W)
#define CELULA_K  4   // Queratinócitos (K)

// Estrutura das células
#define CELULAS_SIZEOF   40
#define CELULAS_SINGLE_CELL_SIZEOF   8
#define CELULAS_L_OFFSET  0
#define CELULAS_M_OFFSET  8
#define CELULAS_S_OFFSET  16
#define CELULAS_W_OFFSET  24
#define CELULAS_K_OFFSET  32
#define CELULAS_POSICAO_OR_OFFSET 0  
#define CELULAS_POSICAO_XP_OFFSET 1  
#define CELULAS_POSICAO_XM_OFFSET 2  
#define CELULAS_POSICAO_YP_OFFSET 3  
#define CELULAS_POSICAO_YM_OFFSET 4  
#define CELULAS_POSICAO_ZP_OFFSET 5  
#define CELULAS_POSICAO_ZM_OFFSET 6  
#define CELULAS_NOVO_VALOR_OFFSET 7

// Informações de acesso à estrutura "parametrosMalhaGPU"
#define OFFSET_COMPUTACAO       0
#define LENGTH_COMPUTACAO       1
#define MALHA_DIMENSAO_X        2
#define MALHA_DIMENSAO_Y        3
#define MALHA_DIMENSAO_Z        4
#define MALHA_DIMENSAO_POSICAO_Z 5
#define MALHA_DIMENSAO_POSICAO_Y 6
#define MALHA_DIMENSAO_POSICAO_X 7
#define MALHA_DIMENSAO_CELULAS  8
#define NUMERO_PARAMETROS_MALHA 9

// Constantes biológicas e fracionárias
__constant float deltaT = 1e-6;
__constant float deltaX = 0.5, deltaY = 0.5, deltaZ = 0.5;
__constant float r1 = 0.5, r2 = 0.4;  // Taxas de crescimento para L e M
__constant float k1 = 5.0, k2 = 4.0;  // Capacidade de suporte para L e M
__constant float gamma1 = 0.065, gamma2 = 0.05, gamma3 = 3.0;
__constant float delta1 = 0.01, delta2 = 0.01, delta3 = 0.01; // Parâmetros específicos de queratinócitos e citocinas
__constant float beta1 = 0.065;
__constant float mu1 = 0.07, mu2 = 0.02, mu3 = 0.7;
__constant float theta1 = 0.5, theta2 = 0.7, theta3 = 0.5;
__constant float lambda = 0.4;
__constant float Vt = 1.0f;
__constant float xi = 0.075; // Parâmetro específico para queratinócitos

// Função Laplaciano
float Laplaciano(int celulaOffset, float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
	return ((xPosicaoGlobal > 0 && xPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET])/(deltaX*deltaX) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)-xPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaX*deltaX)) + (xPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaX*deltaX)))) + ((yPosicaoGlobal > 0 && yPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET])/(deltaY*deltaY) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)-yPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaY*deltaY)) + (yPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaY*deltaY)))) + ((zPosicaoGlobal > 0 && zPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET])/(deltaZ*deltaZ) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)-zPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaZ*deltaZ)) + (zPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaZ*deltaZ))));
}

float Quimiotaxia(int celulaOffset, float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
	return (((((xPosicaoGlobal > 0) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XM_OFFSET]) > 0) ? (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET] / deltaX) : (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaX)) : 0) + ((xPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_X] - 1) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaX) : ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] / deltaX)) : 0))/deltaX) + ((((yPosicaoGlobal > 0) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YM_OFFSET]) > 0) ? (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET] / deltaY) : (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaY)) : 0) + ((yPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_Y] - 1) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaY) : ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] / deltaY)) : 0))/deltaY) + ((((zPosicaoGlobal > 0) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) > 0) ? (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET] / deltaZ) : (-(celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaZ)) : 0) + ((zPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_Z] - 1) ? (((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaZ) : ((celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] / deltaZ)) : 0))/deltaZ));
}







void CalcularPontos(float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
    // Células e variáveis
    float L = celulas[CELULAS_L_OFFSET + CELULAS_POSICAO_OR_OFFSET];
    float M = celulas[CELULAS_M_OFFSET + CELULAS_POSICAO_OR_OFFSET];
    float S = celulas[CELULAS_S_OFFSET + CELULAS_POSICAO_OR_OFFSET];
    float W = celulas[CELULAS_W_OFFSET + CELULAS_POSICAO_OR_OFFSET];
    float K = celulas[CELULAS_K_OFFSET + CELULAS_POSICAO_OR_OFFSET];

    // Quimiotaxia em L e M em direção a S
    float quimiotaxiaL = Quimiotaxia(CELULAS_L_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
    float quimiotaxiaM = Quimiotaxia(CELULAS_M_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);

    // Difusão em L, M, S e W
    float difusaoL = Laplaciano(CELULAS_L_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
    float difusaoM = Laplaciano(CELULAS_M_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
    float difusaoS = Laplaciano(CELULAS_S_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
    float difusaoW = Laplaciano(CELULAS_W_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);

    // Atualização das células L, M, S e W com base nas equações
    celulas[CELULAS_L_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = max(0.0f, ((-mu1 * L) + (L * r1 * (1 - L / k1)) + (gamma1 * L * S) - (delta1 * (1 - W) * L * M) + difusaoL) * deltaT + L);
    
    celulas[CELULAS_M_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = max(0.0f, ((-mu2 * M) + (M * r2 * (1 - M / k2)) - (beta1 * (1 - W) * L * M) + difusaoM) * deltaT + M);
    
    celulas[CELULAS_S_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = max(0.0f, (((theta1 + W) * L * S) / (theta2 + S + theta3 * S) - mu3 * S + difusaoS) * deltaT + S);
    
    celulas[CELULAS_W_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = max(0.0f, ((-gamma3 * W) + Vt) * deltaT + W);
    
    // K (pode precisar ser ajustado conforme a lógica específica do modelo)
    celulas[CELULAS_K_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = max(0.0f, ((xi * (1 - W) * L * M) + (delta3 * (1 - W) * L * K) - lambda * K) * deltaT + K);
}

__kernel void ProcessarPontos(__global float *malhaPrincipalAtual, __global float *malhaPrincipalAnterior, __constant int *parametrosMalhaGPU)
{
    int globalThreadID = get_global_id(0);
    
    float celulas[CELULAS_SIZEOF];

    // Descobrir posição 3D local na malha.
    int posicaoGlobalZ = (globalThreadID / (parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X]));
    int posicaoGlobalY = (globalThreadID % (parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X])) / parametrosMalhaGPU[MALHA_DIMENSAO_X];
    int posicaoGlobalX = (globalThreadID % (parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X])) % parametrosMalhaGPU[MALHA_DIMENSAO_X];

    if (posicaoGlobalZ >= parametrosMalhaGPU[MALHA_DIMENSAO_Z])
    {
        return;
    }

    //**************************************
    // Preencher células para calcular EDOs.
    //**************************************

    int malhaIndex = (posicaoGlobalZ * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z]) + 
                     (posicaoGlobalY * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y]) + 
                     (posicaoGlobalX * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X]);

    // Loop por todas as células.
    for (int count = 0; count < MALHA_TOTAL_CELULAS; count++)
    {   

        // Origem.
        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_OR_OFFSET] = malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS])];

        // Vizinhança.
        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_ZP_OFFSET] = ((posicaoGlobalZ + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_Z])) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z])] : 0.0f;

        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_ZM_OFFSET] = ((posicaoGlobalZ - 1 >= 0)) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (-1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z])] : 0.0f;

        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_YP_OFFSET] = ((posicaoGlobalY + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_Y])) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y])] : 0.0f;

        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_YM_OFFSET] = ((posicaoGlobalY - 1 >= 0)) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (-1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y])] : 0.0f;

        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_XP_OFFSET] = ((posicaoGlobalX + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_X])) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X])] : 0.0f;

        celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_XM_OFFSET] = ((posicaoGlobalX - 1 >= 0)) ? 
            malhaPrincipalAnterior[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + (-1 * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X])] : 0.0f;
    }

    //**************************************
    // Atualizar malha com pontos calculados.
    //**************************************
    
    CalcularPontos(celulas, posicaoGlobalX, posicaoGlobalY, posicaoGlobalZ, parametrosMalhaGPU);

    // Loop por todas as células para atualizar a malha atual.
    for (int count = 0; count < MALHA_TOTAL_CELULAS; count++)
    {   

    

        malhaPrincipalAtual[malhaIndex + ((CELULA_L + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS])] = celulas[((CELULA_L + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_NOVO_VALOR_OFFSET];
    }
}
