#include <cstdlib>
#include <vector>
#include <climits>
#include <math.h>
//#include <matplot/matplot.h>
using point = std::pair<double, double>;
using namespace std;

#include<iostream>
#include <locale.h>
#include <cstring>
#include <ctime>
#include <bits/stdc++.h>
#include <algorithm>

using namespace std;
#define tamMatriz 5


//https://www.tutorialspoint.com/compile_cpp_online.php

class Grafo {
    int tamGrafo;
    vector<vector<int>> verticeAdjacente;

public:

    Grafo(int v) :
            tamGrafo(v), verticeAdjacente(v) {
    }

    void adicionaVertice(int no, int vertice) {
        verticeAdjacente[no].push_back(vertice);
    }







    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------  fim DSF               -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------

    void depth_first(int s, vector<bool> &visited) {
        visited[s] = true;
        cout << s << " ";
        for (int u : verticeAdjacente[s]) {
            if (!visited[u])
                depth_first(u, visited);
        }
    }

//Na teoria dos grafos, busca em profundidade (ou busca em profundidade-primeiro, também conhecido em inglês por
// Depth-First Search - DFS) é um algoritmo usado para realizar uma busca ou travessia numa árvore, estrutura de árvore
// ou grafo. Intuitivamente, o algoritmo começa num nó raiz (selecionando algum nó como sendo o raiz, no caso de um grafo)
// e explora tanto quanto possível cada um dos seus ramos, antes de retroceder(backtracking).
//Thus, for the adjacencymatrix representation, the traversal time is in O(V^2),and for the adjacency listrepresentation,
// it  is  in O(|V|+|E|) where|V|and|E| are  the  number  of  thegraph’s vertices and edges, respectively.
    void DFS() {
        //queue<int>fila;
        vector<bool> visitado(tamGrafo); // guarda o vetor visitado.

        //inicializa o vetor de 0 a N => false
        for (int i = 0; i < tamGrafo; i++) {
            visitado[i] = false;
        }

        //   visitado[no_raiz] = true;
        // cout << no_raiz << " ";



        cout << "DFS: ";
        for (int i = 0; i < tamGrafo; ++i) {
            if (!visitado[i]) {
                depth_first(i, visitado);
            }
        }

    }

    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------  fim BSF               -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------

    /*

    no_atual=raiz
    visitado[no_atual]=true

    enquanto fila != vazio:

         percorrer todos os vizinhos
             enquanto a fila != vazio
                  v=fila.frente();
                  fila.remove
                  se nó não foi visitado
                     fila.insere(v)
                     visitado(v) =true

     */
//Considerando um grafo representado em listas de adjacência, o pior caso, aquele em que todos os vértices e arestas
// são explorados pelo algoritmo,
// Thus, for the adjacencymatrix representation, the traversal time is in O(V^2),and for the adjacency listrepresentation,
// it  is  in O(|V|+|E|) where|V|and|E| are  the  number  of  thegraph’s vertices and edges, respectively.
//  que significa o número de operações sobre todos os vértices que possui uma complexidade constante O ( 1 )
//  para cada vértice uma vez que t.o.d.o vértice é enfileirado e desenfileirado uma unica vez.
    void BFS(int no_raiz) {
        queue<int> fila;
        vector<bool> visitado(tamGrafo);

        for (int i = 0; i < verticeAdjacente.size(); i++) {
            visitado[i] = false;
        }

        fila.push(no_raiz); // insere no final da fila
        visitado[no_raiz] = true;

        cout << "BFS - Visitando: " << no_raiz << endl;
        while (!fila.empty()) {
            int vertice = fila.front(); //pega o primeiro elemento
            fila.pop(); // remove elemento da frente da fila
            cout << vertice << "...";
            for (int no : verticeAdjacente[vertice]) {
                if (!visitado[no]) {
                    fila.push(no);
                    visitado[no] = true;
                }
            }
        }
        cout << endl;
    }

};







//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------  fim caixeiroViajante  -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//https://github.com/adriano-op/paa104_lista_2.git


template<typename It>
bool proximaPermutacao(It begin, It end) {
    if (begin == end)
        return false;

    It i = begin;
    ++i;
    if (i == end)
        return false;

    i = end;
    --i;

    while (true) {
        It j = i;
        --i;

        if (*i < *j) {

            It k = end;

            while (!(*i < *--k))
                /* pass */;

            iter_swap(i, k);

            reverse(j, end);
            return true;
        }

        if (i == begin) {
            reverse(begin, end);
            return false;
        }
    }
}

void troca(vector<int> vetor, int i, int j) {
    int aux = vetor[i];
    vetor[i] = vetor[j];
    vetor[j] = aux;
}

int permuta(vector<int> vetor, int inf, int sup) {
    if (inf == sup) {
        for (int i = 0; i <= sup; i++) {
            vetor.push_back(i);
            std::cout << i << "  ";
        }

        std::cout << "\n" << endl;
    } else {
        for (int i = inf; i <= sup; i++) {
            troca(vetor, inf, i);
            permuta(vetor, inf + 1, sup);
            troca(vetor, inf, i); // backtracking
        }
    }
}
//O problema pode ser convenientemente modelado por um gráfico ponderado, com os vértices do gráfico representando as
// cidades e os pesos das arestas especificando as distâncias. Então, o problema pode ser classificado como o problema
// de encontrar o menor circuito de Hamiltoniano do gráfico

// O(n!)
int caixeiroViajante(int matrizV[tamMatriz][tamMatriz]) {
    vector<int> vizinho; // armazena as permutações da cidades há serem visitadas.
    int aux = 0;
    int pesoMinCaminho = 1000000000; // armazenar peso mínimo do menor caminho.
    int countCaminhoTotal = 0, countMenorCaminho = 0;

    for (int i = 0; i < tamMatriz - 1; i++) {
        vizinho.push_back(i);
    }

    do {
        int somaPesoCaminhoAtual = 0;
        int x = aux;

        for (int i = 0; i < vizinho.size(); i++) {
            cout << "   " << vizinho[i]; // imprime o vetor
        }

        cout << "   " << endl;

        for (int i = 0; i < vizinho.size(); i++) {
            somaPesoCaminhoAtual = 0;
            int y = vizinho[i];

            cout << "[ " << x << " , " << y << " ]  matrizV[x][y]:  " << matrizV[x][y];
            cout << " " << endl;

            somaPesoCaminhoAtual += matrizV[x][y];
            cout << " soma Peso Caminho: " << somaPesoCaminhoAtual << endl;
            x = vizinho[i];
        }
        // somaPesoCaminhoAtual += matrizV[x][aux]; // matriz[x][0]


        if (pesoMinCaminho > somaPesoCaminhoAtual) {
            pesoMinCaminho = somaPesoCaminhoAtual;
            countMenorCaminho++;
        }

        countCaminhoTotal++;
    } while (proximaPermutacao(vizinho.begin(), vizinho.end()));


    cout << "Total de caminhos percorridos: " << countCaminhoTotal;
    cout << "\n O menor peso: " << pesoMinCaminho << " ,no caminho: " << countMenorCaminho << endl;
    return pesoMinCaminho;
}

//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------  fim problemaDaMochila -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------

//O problema da mochila é um problema de otimização combinatória: dado um conjunto de itens, cada um com um peso e
// um valor, determine o número de cada item a incluir em uma coleção de modo que o peso total seja menor ou igual a
// um determinado limite e o valor total é o maior possível.
//O(n).
int knaps(int tamanhoMochila, int peso[], int valor[], int capacidade, int i) {
    int somarValor, contaValor;

    if (i == tamanhoMochila || capacidade <= 0) {
        return 0;
    }


    if (peso[i] < capacidade) { // peso menor que a capacidade
        somarValor = knaps(tamanhoMochila, peso, valor, capacidade - peso[i], i + 1) + valor[i];
        contaValor = knaps(tamanhoMochila, peso, valor, capacidade, i + 1);
        int m = max(somarValor, contaValor);
        return m;
    } else {
        contaValor = knaps(tamanhoMochila, peso, valor, capacidade, i + 1);
        return contaValor;
    }
}


//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------  fim convexHull        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------

struct Point {
    int x, y;
};

int orientacao(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;
    return (val > 0) ? 1 : 2; // no sentido horário ou anti-horário
}


//Um conjunto de pontos (finitos ou infinitos) no plano é chamado de convexo se, para quaisquer dois pontos p e q no
// conjunto, t.o.d.o o segmento de linha com as extremidades em p e q pertence ao conjunto.
//What  is  the  time  efficiency  of  this  algorithm?  It  is  in O(n^3)

void convexHull(Point points[], int n) {
    // Deve haver pelo menos 3 pontos
    if (n < 3) {
        cout << " Deve haver pelo menos 3 pontos " << endl;
        return;
    }

    vector<Point> hull; //

    // Encontra o ponto mais à esquerda
    int left = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].x < points[left].x) {
            left = i;
        }
    }

    // Comece do ponto mais à esquerda, continue movendo no sentido anti-horário
    // até chegar ao ponto inicial novamente. Este loop é executado O (h)
    // vezes em que h é o número de pontos no resultado ou na saída.
    int p = left, q;
    do {
        // Adiciona o ponto atual ao resultado
        hull.push_back(points[p]);

        // Procure um ponto 'q' tal que a orientação (p, q,
        // x) é anti-horário para todos os pontos 'x'.
        // A ideia é manter um registro do último ponto mais visitado no sentido anti-horário em q.
        //Se qualquer ponto 'i' for mais anti-horário do que q, atualize q.
        q = (p + 1) % n;
        for (int i = 0; i < n; i++) {
            // // Se i for mais anti-horário do que o q atual, então // atualize o q
            if (orientacao(points[p], points[i], points[q]) == 2) {
                q = i;
            }
        }

        // // Agora q é o mais anti-horário em relação a p, defina p como q para a próxima iteração,
        // de modo que q seja adicionado ao resultado 'hull'
        p = q;

    } while (p != left);

    for (int i = 0; i < hull.size(); i++) {
        cout << "(" << hull[i].x << ", " << hull[i].y << ")\n";
    }
}

//  dada uma string de n caracteres chamada
// de texto e uma string de m caracteres (m ≤ n) chamada de padrão, encontre uma substring do texto que corresponda ao padrão.
// Para ser mais preciso, queremos findi - o índice do caractere mais à esquerda da primeira substring correspondente

// O(n)
int bruteForceStringMatch(char palavra[], char termo[]) {
    // baseado no ALGORITHM BruteForceStringMatch(T[0..n−1],P[0..m−1]
    cout << " palavra: " << strlen(palavra) << endl;
    cout << " termo : " << strlen(termo) << endl;

    if (strlen(palavra) < strlen(termo))
        return -1;

    if (strlen(palavra) == 0 || strlen(termo) == 0) {
        return -1;
    }
    for (int i = 0; i < strlen(palavra) - strlen(termo) + 1; i++) {
        int j = 0;
        while (j < strlen(termo) && palavra[i + j] == termo[j]) {
            j++;
        }
        if (j == strlen(termo)) {
            return i;
        }

    }
    return -1;
}

//------------------------                            -------------------------------
//------------------------                            -------------------------------
//------------------------  fim Imprimir Vetor        -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------

void ImprimeVetor(int v[], int Tam) {

    for (int i = 0; i < Tam; i++) {
        cout << v[i] << " | ";
    }

    cout << endl;
}
// ------------------------                            -------------------------------

// ------------------------                            -------------------------------
// ------------------------  Fim bruteForceClosestPair -------------------------------
// ------------------------                            -------------------------------
// ------------------------                            -------------------------------

double bruteForceClosestPair(vector<point> vetorPair) {
    //O algoritmo abaixo calcula a distância entre os dois pontos mais próximos;
    //Baseado no ALGORITHM BruteForceClosestPair(P)

    cout << "Tamnho vetor= " << vetorPair.size() << endl;

    if (vetorPair.size() < 2) {
        cout << "O vetor deve possuir no minimo dois elemntos." << endl;
        return 0;
    }

    for (int i = 0; i < vetorPair.size() - 1; i++) {
        cout << "points [ " << vetorPair[i].first << ", " << vetorPair[i].second << "]" << endl;
    }


    double distancia_min = __DBL_MAX__; //  é o número máximo de ponto flutuante finito representável.

    for (int i = 0; i < vetorPair.size() - 1; i++) {
        for (int j = i + 1; j < vetorPair.size(); j++) {
            //calcula a distancia min da subtração do quadrado entre dois pontos.
            double p1 = pow(vetorPair[i].first - vetorPair[j].first, 2) +
                        pow(vetorPair[i].second - vetorPair[j].second, 2);

            if (distancia_min > p1) {
                distancia_min = p1;
            }
        }
    }

    return sqrt(distancia_min); //calcula a raiz
}


//------------------------                            -------------------------------
//------------------------                            -------------------------------
//------------------------  fim SequentialSearch2     -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------

//O(n)
// Para repetir, o algoritmo simplesmente compara elementos sucessivos de uma determinada lista com uma determinada
// chave de pesquisa até que uma correspondência seja encontrada (pesquisa bem-sucedida) ou a lista se esgote sem
// encontrar uma correspondência (pesquisa malsucedida).
int SequentialSearch2(int v[], int k, int Tam) {
    int i = 0;

    while (v[i] != k) {
        i++;
    }
    if (i <= Tam) {
        return i;
    } else {
        return -1;
    }
}

//Thus, selection sort is a O(n^2)algorithm on all inputs.
//Começamos a ordenação por seleção examinando toda a lista fornecida para encontrar seu menor elemento e trocá-lo pelo
// primeiro elemento, colocando o menor elemento em sua posição final na lista ordenada.
// Em seguida, examinamos a lista, começando com o segundo elemento, para encontrar o menor entre os últimos n − 1
// elementos e trocá-lo pelo segundo elemento, colocando o segundo menor elemento em sua posição final
void selectionSort(int v[], int n) {
    int i, j, min, aux;
    for (i = 0; i < n - 2; i++) {
        min = i;
        for (j = i + 1; j < n - 1; j++)
            if (v[j] < v[min]) {
                min = j;
            }

        aux = v[min];
        v[min] = v[i];
        v[i] = aux;
    }

    cout << "\n SelectionSort: Elementos do array em ordem crescente:\n" << endl;
    ImprimeVetor(v, n);
}
// ------------------------  Inicio selectionSort -------------------------------


//Another brute-force application to the sorting problem is to compare adjacentelements  of  the  list  and  exchange
// them  if  they  are  out  of  order.
// O(n^2).
void Bubblesort(int v[], int TAM) {
    //baseado no ALGORITHM BubbleSort(A[0..n−1]),
    //no livro /Introduction to the Design and Analysis of Algorithms (3rd ed.) [Levitin 2011-10-09].pdf

    for (int j = 0; j < TAM - 2; j++) {
        for (int i = 0; i < TAM - j - 2; i++) {

            if (v[i + 1] < v[i]) {
                int aux = v[i];
                v[i] = v[i + 1];
                v[i + 1] = aux;
            }
        }
    }


    cout << "\n Bubblesort: Elementos do array em ordem crescente:\n" << endl;
    ImprimeVetor(v, TAM);
    cout << " " << endl;
}

void inicializaVetor(int v[], int tam) {
    unsigned seed = time(0);
    srand(seed);

    for (int i = 0; i < tam; i++) {
        v[i] = 1 + rand() % 1000;
        //     cout << v[i] << " ";
    }
}

void gp(vector<int> current_path, vector<int> to_do) {
    if (to_do.empty()) {
        //       inicializaVetor(current_path, current_path.size());
    } else {
        for (int i = 0; i < to_do.size(); i++) {
            vector<int> path(current_path);
            path.push_back(to_do[i]);
            vector<int> left_to_do(to_do);
            left_to_do.erase(left_to_do.begin() + i);

            gp(path, left_to_do);
        }
    }
}

void generation(int n) {
    vector<int> current_path = {0}; // tirando o ), ele gera para todas as cidades
    vector<int> to_do;

    for (int i = 0; i < n; i++) {
        to_do.push_back(i);
    }

    gp(current_path, to_do);
}

void print_vector(std::vector<int> &v) {
    for (auto element : v) {
        std::cout << element << ' ';
    }

    std::cout << std::endl;
}

void print_segments(std::vector<std::vector<int>> &segments) {
    for (auto segment : segments) {
        std::cout << "[" << segment[0] << " , " << segment[1] << "] " << " [" << segment[2] << " , " << segment[3]
                  << "]" << std::endl;
    }
}

std::vector<int> get_vals(int n) {
    std::vector<int> v;
    for (size_t i = 0; i < n; i++) {
        v.push_back(rand() % 1000);
    }

    return v;
}


//void plot_stuff(std::vector<int>& x, std::vector<int>& y , std::vector<std::vector<int>>& segments){
//    auto fig = matplot::plot(x, y, "o");
//    fig->marker_size(10);
//    fig->marker_face_color("#0000FF");
//
//    matplot::hold(matplot::on);
//    for (auto s : segments) {
//        matplot::plot(std::vector<int>{s[0],s[2]}, std::vector<int>{s[1],s[3] }, "-")->line_width(2);
//    }
//
//    matplot::show();
//
//}

std::vector<std::vector<int>> get_convex_hull(std::vector<int> &x, std::vector<int> &y) {
    //On³ complexidade

    std::vector<vector<int>> segments;
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = i + 1; j < x.size(); j++) {
            int a = y[j] - y[i];
            int b = x[i] - x[j];
            int c = x[i] * y[j] - y[i] * x[j];

            // se ax +by+c < 0 = soma -1
            // se não soma +1
            int sum_sign = 0;
            for (size_t k = 0; k < x.size(); k++) {
                if (k != i && k != j) {
                    //ax + by -c = 0
                    if ((a * x[k] + b * y[k] - c) < 0) {
                        sum_sign--;
                    } else {
                        sum_sign++;
                    }
                }
                if (std::abs(sum_sign) == (x.size() - 2)) {
                    segments.push_back({x[i], y[i], x[j], y[j]});
                }
            }
        }

    }
    return segments;

}


int main(int argc, char **argv) {
//1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
//1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
//------------------------------------  1 > Inicio  Bubblesort            ----------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------


//    cout << "Algoritmo N° 1: Bubblesort" << endl;
//    cout << "Insira um valor, para que função gere as entradas do sistema: ";
//    int TAM;
//    cin >> TAM;
//
//    int buble[TAM];
//    inicializaVetor(buble, TAM);
//
//    cout << "array size: " << TAM << endl;
//    cout << "Função Bubblesort \n";
//    cout << "\n Ordem atual dos itens no array:\n" << endl;
//    ImprimeVetor(buble, TAM);
//    Bubblesort(buble, TAM);






//2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
//2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//---------------------------------  2 > Inicio  selectionSort         -------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------



//    cout << "Algoritmo N° 2: selectionSort" << endl;
//    cout << "Insira um valor, para que função gere as entradas do sistema: ";
//    int TAM;
//    cin >> TAM;
//
//    if (TAM == 0) {
//        cout << "\n N==0,  programa será encerrado. " << endl;
//        return 0;
//    }
//
//    int select[TAM];
//    inicializaVetor(select, TAM);
//    cout << "array size: " << TAM << "\n";
//    cout << "Função selectionSort \n";
//    cout << "Ordem atual dos itens no array:\n" << endl;
//    ImprimeVetor(select, TAM);
//    selectionSort(select, TAM);




//3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
//----------------------------------------------------------------------------------------------------------------------
//-----------------------------  3 > Inicio  SequentialSearch2     -----------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------

//    cout << "Algoritmo N° 3: SequentialSearch2" << endl;
//    cout << "Insira um valor, para que função gere as entradas do sistema: ";
//    int TAM;
//    cin >> TAM;
//
//    if (TAM == 0) {
//        cout << "\n N==0,  programa será encerrado. " << endl;
//        return 0;
//    }
//
//    int sequential[TAM];
//    inicializaVetor(sequential, TAM);
//    cout << "Função SequentialSearch2 \n" << endl;
//    cout << "Ordem atual dos itens no array:\n" << endl;
//
//    ImprimeVetor(sequential, TAM);
//
//    cout << "\n Informe o elemento para a busca: ";
//    int elemento;
//    cin >> elemento;
//    int w = SequentialSearch2(sequential, elemento, TAM);
//
//    if (w != -1) {
//        cout << "Elemento " << elemento << " encontrado na posição: " << w << endl;
//    } else {
//        cout << "Elemento não encontrado " << endl;
//    }




//4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
//----------------------------------------------------------------------------------------------------------------------
//--------------------------------  4 > Inicio  bruteForceStringMatch --------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------



//    cout << "Algoritmo N° 4: bruteForceStringMatch" << endl;
//
//    char letras[50];
//    char palavra[50];
//
//    cout << "Insira uma palavra: " << endl;
//    cin.getline(palavra, 50);
//
//    cout << " Insira o termo da busca : " << endl;
//    cin.getline(termo, 50);
//
//    int index = bruteForceStringMatch(palavra, termo);
//
//
//    if (index == -1) {
//        cout << "Termo não encontrado" << endl;
//    }
//    else {
//        cout << "Termo encontrado na posição: " << index << endl;
//    }




//5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
//5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
//------------------------                            -------------------------------
//------------------------  5 > bruteForceClosestPair -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------
//------------------------                            -------------------------------

// Ao iniciar a execução da função, o sistema gera o tamanho do vetor e inicializa automaticamente
// com as coordenadas aleatorias de 0 a 100.
// de posse do vector de coordenadas, a função calcula a distância entre os dois pontos mais próximos;


//    cout << "Algoritmo N°5: BruteForceClosestPair" << endl;
//    unsigned seed = time(0);
//    srand(seed);
//
//    int tam = rand() % 100;
//    vector<point> pair;
//    for (int i = 0; i < tam; i++)
//        pair.push_back(point((double) (rand() % 100), (double) (rand() % 100)));
//
//    cout << bruteForceClosestPair(pair) << endl;



//6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
//6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//------------------------  6 > convexHull        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------

//        cout << "Algoritmo N° 6: convexHull" << endl;

// O struct points fornce um conjunto de pontos, que serão passados para o convex_hull calcular o segmento de linha com
// os pontos finais em p e q pertence ao conjunto



//    Point points[] = {
//            {0, 3},
//            {2, 2},
//            {1, 1},
//            {2, 1},
//            {3, 0},
//            {0, 0},
//            {3, 3}
//    };
//
//    int tam = sizeof(points) / sizeof(points[0]);
//    convexHull(points, tam);


    //   https://github.com/adriano-op/paa104_lista_2_adriano_andrade.git

    // exemplo do professor

//    srand((unsigned int) time(0));
//    int n = 15;
//    std::vector<int> x = get_vals(n);
//    std::vector<int> y = get_vals(n);
////
//    std::vector<std::vector<int>> segments= get_convex_hull(x,y);
//    print_segments(segments);
//    //   plot_suttfy();
//
//    generation(4);




//7777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
//7777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
//----------------------------------------------------------------------------------------------------------------------
//---------------------------------------  7 > Caixeiro Viajante -------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------


//    int matriz[tamMatriz][tamMatriz] = {
//            {0,  2,  4,  6,  8},
//            {10, 0,  14, 16, 18},
//            {20, 22, 0,  26, 28},
//            {30, 32, 34, 0,  38},
//            {40, 42, 44, 46, 0}
//    };
//
//
//    caixeiroViajante(matriz);



//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------  8 > problemaDaMochila -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------
    //------------------------                        -------------------------------

//    cout << "Algoritmo N° 8: problemaDaMochila" << endl;
//
//    int n = 5;
//    int peso[5] = {2, 4, 3, 5, 5};
//    int valor[5] = {3, 4, 1, 2, 6};
//    int capacidade = 12;
//
//    cout << "Capacidade: " << capacidade << endl;
//    cout << "Valor maximo da mochila: " << knaps(n, peso, valor, capacidade, 0) << endl;


//9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
//9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
//------------------------                         -------------------------------
//------------------------                        -------------------------------
//------------------------  9 >      BFS          -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------


//                     0
//                /            \
    //              1               2
//          /    \          /     \
    //        3       4      5          6
//                                   \
    //                                    7


//    Grafo grafo(8);
//    grafo.adicionaVertice(0, 1);
//    grafo.adicionaVertice(0, 2);
//    grafo.adicionaVertice(1, 3);
//    grafo.adicionaVertice(1, 4);
//    grafo.adicionaVertice(2, 5);
//    grafo.adicionaVertice(2, 6);
//    grafo.adicionaVertice(6, 7);
//


//    grafo.BFS(0); // nó raiz = 0






//------------------------                         -------------------------------
//------------------------                        -------------------------------
//------------------------  10 >      DFS          -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------
//------------------------                        -------------------------------

//                     0
//                /            \
    //              1               2
//          /    \          /     \
    //        3       4      5          6
//                                   \
    //                                    7
/*
    Grafo grafinho(8);
    grafinho.adicionaVertice(0, 1);
    grafinho.adicionaVertice(0, 2);
    grafinho.adicionaVertice(1, 3);
    grafinho.adicionaVertice(1, 4);
    grafinho.adicionaVertice(2, 5);
    grafinho.adicionaVertice(2, 6);
    grafinho.adicionaVertice(6, 7);



    grafinho.DFS(); // nó raiz = 0

 */
    return 0;
}


