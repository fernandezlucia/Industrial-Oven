#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<string>
#include<chrono>

using namespace std;

double radioInterno;
double radioExterno;
double cantRadios;
double cantAngulos;
vector<double> radios;
vector<double> angulos;
float iso;
int ninst;
vector<double> t_Int;
vector<double> t_Ext;
int cantFilas;
int cantColumnas;
double stepRadio;
double stepAngulo;
vector<vector<double>> instancias;

vector<double> resolucionSuperior(vector<vector<double>> triangulada){ 
    vector<double> x(triangulada.size());
    int columna_b = triangulada[0].size()-1;
    for(int i = triangulada.size()-1; i >= 0; --i){
        double sum = 0;
        for(int j = i+1; j < columna_b; ++j)
            sum += triangulada[i][j]*x[j]; 
        x[i] = (triangulada[i][columna_b] - sum)/triangulada[i][i];
    }
    return x;
}

vector<double> resolucionInferior(vector<vector<double>> triangulada){
    vector<double> x(triangulada.size());
    int columna_b = triangulada[0].size()-1;

    for(int i = 0; i <= triangulada.size()-1; i++){
        double sum = 0;
        for(int j = 0; j <= i; j++)
            sum += triangulada[i][j]*x[j]; 
        x[i] = (triangulada[i][columna_b] - sum)/triangulada[i][i];
    }
    return x;
}

vector<vector<double>> L(vector<vector<double>> matrix){
vector<vector<double>> res(cantFilas, vector<double>(cantColumnas));
//Para cada columna, una vez que guardamos el coeficiente tenemos que restar la fila de matrix  y esto se repite en todos los elementos debajo de la diagonal
    for (int i = 0; i < cantFilas; i++) {
        for (int j = i + 1; j < cantFilas; j++){ 
            double coef = matrix[j][i] / matrix[i][i];
            res[j][i] = coef;
            for (int k = 0; k < cantFilas ; k++){
             if(coef != 0) 
                matrix[j][k] -= coef * matrix[i][k];
            }
        }
    }
    for(int i=0; i< cantFilas; i++){
        res[i][i] = 1;
    }
    return res;
}

vector<double> resoLU(vector<vector<double>> L, vector<vector<double>> U, vector<double> b){
    vector<double> y;
    vector<double> res;
    for(int i=0;i<L.size();i++)
        L[i].push_back(b[i]);
    y = resolucionInferior(L);
    for(int i=0;i<U.size();i++)
        U[i].push_back(y[i]);
    res = resolucionSuperior(U);
    return res;
}

void procesarArchivo(string nombreArchivo){
    ifstream entrada;
    entrada.open(nombreArchivo);

    entrada >> radioInterno >> radioExterno >> cantRadios >> cantAngulos >> iso >> ninst;

    vector<vector<double>> resInst(ninst, vector<double>(2*cantAngulos));
    instancias = resInst;
  

    for(int i = 0; i < ninst; i++){
        for(int j = 0; j < 2*cantAngulos; j++){
            entrada >> instancias[i][j];
        }
    }
   

    t_Int = vector<double>(cantAngulos);
    t_Ext = vector<double>(cantAngulos);
    angulos = vector<double>(cantAngulos);
    radios = vector<double>(cantRadios);

    entrada.close();

    cantFilas = (cantRadios) * cantAngulos;
    cantColumnas = (cantRadios) * cantAngulos;

    stepRadio = (radioExterno - radioInterno)/(cantRadios-1);
    stepAngulo = 2/cantAngulos*M_PI;

    for(int i = 0; i < cantAngulos; i++){
        angulos[i] = stepAngulo*i;
    }
     for(int i = 0; i < cantRadios; i++){
        radios[i] = radioInterno + i*stepRadio;
    }
}


void setearInstancia(int n){
    for(int i = 0; i < cantAngulos; i++)
        t_Int[i] = instancias[n][i];
    

    for(int i = cantAngulos; i < cantAngulos*2; i++)
        t_Ext[i-cantAngulos] = instancias[n][i];
}


vector<double> haceCuentita(int k){
    vector<double> respuesta;
    double deltaR = stepRadio;
    double deltaTita = stepAngulo;

        respuesta.push_back((1/pow(deltaR, 2))-(1/(radios[k % int(cantRadios)]*deltaR)));                                              // j-1,k
        respuesta.push_back(1/(pow(radios[k % int(cantRadios)],2)*pow(deltaTita, 2)));                                                   // j,k+1
        respuesta.push_back((1/(radios[k % int(cantRadios)]*deltaR)) - 2/pow(deltaR, 2) - 2/(pow(radios[k % int(cantRadios)], 2)*pow(deltaTita, 2)));
        respuesta.push_back(1/(pow(radios[k % int(cantRadios)],2)*pow(deltaTita, 2)));                                                 // j, k-1
        respuesta.push_back(1/(pow(deltaR, 2)));                                                                     // j+1,k

    return respuesta;

}

vector<vector<double>> crearMatriz(){     
    vector<vector<double>> res(cantFilas, vector<double>(cantColumnas));

    for(int i = 0; i < cantFilas; i++){      
                                                 //haceCuentita... 0=j-1,k ; 1=j,k+1 ; 2=j,k ; 3=j, k-1 ; 4=j+1, k
       if(i < cantAngulos)
            res[i][i] = 1;

        else if(i > cantFilas-cantAngulos-1)
            res[i][i] = 1;

        else{ 
            //caso general
            vector<double> cuentita = haceCuentita(i);
            res[i][i-cantAngulos] = cuentita[0]; 
            res[i][i-1]  = cuentita[3]; 
            res[i][i] = cuentita[2]; 
            res[i][i+1] = cuentita[1];
            res[i][i+cantAngulos] = cuentita[4];

        }
    }
   return res;
}

vector<vector<double>> extendidaConB (vector<vector<double>> matriz, vector<double> b){
     for (int j = 0; j < b.size(); j++)
        matriz[j].push_back(b[j]);
        
    return matriz;
}

vector<double> vectorB (vector<double> t_int, vector<double> t_ext){
   vector<double> b(cantFilas);
    for(int i = 0; i < cantFilas; i++){
        if(i < cantAngulos)
            b[i] = t_int[i];
        else if(i > cantFilas-cantAngulos-1)
            b[i] = t_ext[i % int(cantAngulos)];
        else
            b[i] = 0;
    }
    return b;
}

vector<vector<double>> eliminacionGaussiana(vector<vector<double>> matrix, vector<double> b){
    double tol = 1e-10;

    for(int i = 0; i < b.size(); i++)
        matrix[i].push_back(b[i]);

    for (int i = 0; i < cantFilas; i++) {
        for (int j = i + 1; j < cantFilas; j++) {
            double coef = matrix[j][i] / matrix[i][i];
            for (int k = i; k < cantColumnas+1; k++) {// el +1 para incluir la columna b
                matrix[j][k] -= coef * matrix[i][k];

                if(abs(matrix[j][k]) < tol)
                    matrix[j][k] = 0;
            }
        }
    }
    return matrix;
}

vector<vector<double>> U(vector<vector<double>> matrix){
    double tol = 1e-10;

    for (int i = 0; i < cantFilas; i++) {
        for (int j = i + 1; j < cantFilas; j++) {
            double coef = matrix[j][i] / matrix[i][i];
            for (int k = i; k < cantColumnas; k++) {
                matrix[j][k] -= coef * matrix[i][k];

                if(abs(matrix[j][k]) < tol)
                    matrix[j][k] = 0;
            }
        }
    }
    return matrix;
}



vector<double> buscoIso(double iso, vector<double> temps){
    vector<double> isoterma(cantAngulos,-1);
    
    for(int i = 0; i < cantAngulos; i++){
        double prev_temp = 0;
        double current_temp = 0;
        for(int j = i; j < temps.size(); j += cantAngulos){
            prev_temp = current_temp;
            current_temp = j;
            if(temps[current_temp] < iso)
                break;
        }
        isoterma[i] = (int)((((temps[prev_temp] - iso) < (iso - temps[current_temp])) ? prev_temp : current_temp) / cantAngulos);
    }
    return isoterma; // posicion = angulo, contenido = radio;
}


void imprimirTemperaturas(vector<double> res, string salida_path){
    ofstream salida;
    
    salida.open(salida_path, ios::app);
        for(int i = 0; i < res.size(); i++){
            salida << res[i] << endl;
                    
        }
    salida.close();
}

void imprimirIsoterma(vector<double> iso, string salida_path){
    ofstream salida;
    
    salida.open(salida_path, ios::app);
        for(int i = 0; i < iso.size(); i++){
            salida << iso[i] << endl;
          
        }
    salida.close();
}



int main(int argc, char **argv){
    if (argc != 4 || atoi(argv[3]) > 1){
        cout << "\nUso incorrecto.\n" << endl;
        cout << "Ejemplo: ./tp test1.in test1.out <metodo>" << endl;
        cout << "Siendo <metodo> uno de los siguientes: " << endl;
        cout << "\t 0: EG." << endl;
        cout << "\t 1: LU." << endl;
        return 1;
    }
    // ejemplo:  $ ./tp1 test1.in test1.out 0
    string entrada_path = argv[1];
    string salida_path = argv[2];
    string temp_path = "temp"; 
    string iso_path = "iso";
    string time_path = "tiempo";
    string lu_path = "LU";
    int metodo = atoi(argv[3]);
    vector<double> res;
    vector<double> isoterma;
    vector<vector<double>> matrizLaplace;
    vector<vector<double>> matrizGauss;
    vector<vector<double>> l;
    vector<vector<double>> u;
    vector<double> b;
    
    switch(metodo){
      case 0: 
      cout << "\n Eliminacion Gaussiana\n" << endl;

            procesarArchivo(entrada_path);
            matrizLaplace = crearMatriz();
            matrizGauss = matrizLaplace;

        for(int j = 0; j < ninst; j++){ 
            setearInstancia(j);
            b = vectorB(t_Int, t_Ext);
            matrizGauss = eliminacionGaussiana(matrizLaplace, b);
            res = resolucionSuperior(matrizGauss); 
            isoterma = buscoIso(iso,res);
            imprimirTemperaturas(res, temp_path + salida_path); 
            imprimirIsoterma(isoterma, iso_path + salida_path);
        }    
        break;
        case 1:
            cout << "\nFactorizacion LU\n" << endl;
            procesarArchivo(entrada_path);
            matrizLaplace = crearMatriz();
            l = L(matrizLaplace);
            u = U(matrizLaplace); 
            for(int j=0; j < ninst; j++){
                setearInstancia(j);
                b = vectorB(t_Int, t_Ext);
                res = resoLU(l,u,b); 
                isoterma = buscoIso(iso,res);
                imprimirTemperaturas(res,lu_path + temp_path + salida_path);
                imprimirIsoterma(isoterma,lu_path + iso_path + salida_path);
                
            } 
            break;
        default:
        cout << "Método inválido" << endl;
        break;
    }
}
