/*
 * ****cparam****
 * (Calculador de parámetros estructurales de nanotubos tipo armchair.)
 * 
 * 
 * Autores:
 * 	Gustavo Dominguez Rodriguez[1]
 * 	Gabriel Ivan Canto Santana [1]
 * 	Jorge Alejandro Tapia Gonzalez [2]
 * 	Cesar Alberto Cab Cauich[2]
 * 
 * 	[1]	Centro de Investigacion en Corrosion, Universidad Autonoma de
 * 		Campeche.
 * 	[2]	Facultad de Ingenieria, Universidad Autonoma de Yucatan
 * 
 * 
 * 
 * ****Cálculo de los parametros estructurales de nanotubos de carbono
 * armchair.****
 * 
 * EL paquete lee el conjunto de bases vectoriales de la celda unitaria,
 * y las posiciones atómicas, en caso de que dichas posiciones sean
 * relativas, las posiciones se calculan se calculan como la combinación
 * lineal de las bases vectoriales, de la forma.
 *		x_i^(c)=X_(i,1) a_1+X_(i,2) a_2+X_(i,3) a_3
 * donde x_i^(c) es la posición en coordenadas cartesianas, que definen
 * la posición del átomo i;  X_(i,1), X_(i,2) y X_(i,3) son los
 * componentes de su vector de coordenadas relativas X_i; mientras que
 * a_1, a_2 y a_3 son los vectores que definen la base de la celda
 * unitaria.
 * 
 * Una vez realizado el procedimiento para cada átomo se procede a
 * calcular el centro de la estructura como la posición media de ellas 
 * [x^(m)], la cual se reporta al usuario, y posteriormente se procede a
 * substraerlo de todas las posiciones atómicas para obtener las
 * posiciones centradas para el átomo i x_i. Ello permitirá simplificar
 * los subsecuentes cálculos.
 * 
 * Se asume que 3 es la dirección axial, y las posiciones atómicas se
 * reordenan de acuerdo a su ángulo en el plano 12 [ϕ_(12,i)], desde 0 a
 * 2π. De la forma
 * 		ϕ_(12,i)=ArcTan[(proy_(x_i) a_2)/(proy_(x_i) a_1 )]
 * donde proy_(x_i) a_1 y proy_(x_i) a_2 son las proyecciones del vector
 * x_i en las direcciones 1 y 2 respectivamente.
 * 
 * Una vez ordenados se empieza a buscar patrones de la forma, dos
 * átomos subsecuentes en posición baja, respecto a 1 y dos átomos en
 * posición alta respecto al mismo. En caso de que se tengan dos átomos
 * en posición alta primero, se procede a rotar los índices
 * identificando a los átomos dos posiciones, para así conseguir que los
 * dos átomos que anteriormente eran los últimos pasen a ser el primero
 * y el segundo y así tener el patrón deseado de dos en posición baja y
 * dos en posición alta. En caso de que se tenga un átomo en posición
 * alta y posteriormente un átomo en posición baja se hace un
 * corrimiento rotativo de índices una posición hacia atrás, para que el
 * que anteriormente era el primer átomo pase a ser el último y
 * nuevamente se tenga el patrón deseado. Y finalmente, si se tiene un
 * átomo en posición baja y un átomo en posición alta, se realiza
 * igualmente un corrimiento rotativo de índices, pero en dirección
 * ascendente, para así tener el átomo que originalmente era el último
 * en la primera posición, y nuevamente, satisfacer el patrón buscado.
 * 
 * Una vez se han centrado las coordenadas y ordenado los índices de 
 * posiciones se procede a calcular los diferentes parámetros a
 * reportar, y después se calcularían las posiciones promedio:
 * 
 * 		Radio. Se calcula simplemente obteniendo la norma de cada
 * posición atómica.
 * 			R_i=‖x_i‖
 * 		Longitud de enlace horizontal. Se calcula como la distancia
 * entre cada átomo con índice impar y el subsecuente. Los índices
 * impares pueden simplemente definirse como 2i-1.
 * 			r^(a)_(2i-1)=‖x_(2i-1)-x_2i‖
 * 
 * 		Longitud de enlace oblicuo. Se calcula como la distancia entre
 * cada átomo con índice par y el subsecuente. Los índices pares pueden
 * simplemente definirse como 2i.
 * 			r^(b)_2i=‖x_2i-x_(2i+1)‖
 * No obstante, también se requiere calcular la distancia entre el átomo
 * par y la posición el átomo siguiente pero la celda unitaria inferior
 * [x↓_(2i+1)].
 * 			r↓^(b)_2i=‖x_2i-x↓_(2i+1)‖
 * En principio, ambas posiciones deben ser idénticas, a no ser que la
 * distancia entre ambos átomos en la dirección 1 no sean exactamente un
 * medio del alto de la celda.
 * 
 * 		Angulo entre enlace horizontal y oblicuo. De este ángulo se
 * pueden calcular dos instancias centradas en cada ángulo. Por lo que
 * habría ocho ángulos por quiralidad.
 * Para las posiciones 1, 5, 9, etc., se calcula el ángulo formado
 * por el átomo del índice anterior (tanto el de la celda central como
 * el de la celda inferior respecto a la dirección 3), el átomo actual y
 * el átomo siguiente.
 * 		θ^(ab)_(4i-3)=ArcCos[({x_(4i-3)-x_(4i-4)}∙{x_(4i-3)-x_(4i-2)})
 * 			/{‖x_(4i-3)-x_(4i-4)‖‖x_(4i-3)-x_(4i-2)‖}],
 * 		θ↓^(ab)_(4i-3)=ArcCos[({x_(4i-3)-x↓_(4i-4)}∙{x_(4i-3)-x_(4i-2)})
 * 			/{‖x_(4i-3)-x↓_(4i-4)‖‖x_(4i-3)-x_(4i-2)‖}]
 * Por otro lado para los ángulos centrados en los átomos 2, 6, 10,
 * etc., se calcula el ángulo formado entre el átomo anterior, el átomo
 * actual y el átomo siguiente (tanto para la celda central como la
 * celda anterior respecto a la dirección 3).
 * 		θ^(ab)_(4i-2)=ArcCos[({x_(4i-2)-x_(4i-3)}∙{x_(4i-2)-x_(4i-1)})
 * 			/{‖x_(4i-2)-x_(4i-3)‖‖x_(4i-2)-x_(4i-1)‖}],
 * 		θ↓^(ab)_(4i-2)=ArcCos[({x_(4i-2)-x_(4i-3)}∙{x_(4i-2)-x↓_(4i-1)})
 * 			/{‖x_(4i-2)-x_(4i-3)‖‖x_(4i-2)-x↓_(4i-1)‖}]
 * En cuanto a los ángulos para los átomos 3, 7, 11, etc., se calcula el
 * ángulo formado entre el átomo anterior, el átomo actual (tanto para
 * la celda central como la celda anterior respecto a la dirección 3) y
 * el átomo siguiente.
 * 		θ^(ab)_(4i-1)=ArcCos[({x_(4i-1)-x_(4i-2)}∙{x_(4i-1)-x_4i})
 * 			/{‖x_(4i-1)-x_(4i-2)‖‖x_(4i-1)-x_4i‖}],
 * 		θ↓^(ab)_(4i-1)=ArcCos[({x↓_(4i-1)-x_(4i-2)}∙{x↓_(4i-1)-x_4i})
 * 			/{‖x↓_(4i-1)-x_(4i-2)‖‖x↓_(4i-1)-x_4i‖}]
 * Finalmente, los ángulos para los átomos 4, 8, 12, etc., se calcula el
 * ángulo formado entre el átomo anterior, el átomo actual y el átomo
 * siguiente (tanto para la celda central como la celda superior
 * respecto a la dirección 3).
 * 		θ^(ab)_4i=ArcCos[({x_4i-x_(4i-1)}∙{x_4i-x_(4i+1)})
 * 			/{‖x_4i-x_(4i-1)‖‖x_4i-x_(4i+1)‖}],
 * 		θ↑^(ab)_4i=ArcCos[({x_4i-x_(4i-1)}∙{x_4i-x↑_(4i+1)})
 * 			/{‖x_4i-x_(4i-1)‖‖x_4i-x↑_(4i+1)‖}]
 * 
 * 		Angulo entre dos enlaces oblicuos. De este ángulo se calcula una
 * sola instancia por cada ángulo. Por lo que habría cuatro ángulos por
 * quiralidad.
 * Para las posiciones 1, 5, 9, etc., se calcula el ángulo formado por
 * el átomo del índice anterior, el átomo actual y el nuevamente átomo
 * del índice anterior pero de la celda inferior respecto a la dirección
 * 3.
 * 		θ^(bb)_(4i-3)=ArcCos[({x_(4i-3)-x_(4i-4)}∙{x_(4i-3)-x↓_(4i-4)})
 * 			/{‖x_(4i-3)-x_(4i-4)‖‖x_(4i-3)-x↓_(4i-4)‖}]
 * Mientras que para las posiciones 2, 6, 10, etc., se calcula el
 * ángulo formado por el átomo del índice siguiente, el átomo actual y
 * el nuevamente átomo siguiente pero de la celda inferior respecto a la
 * dirección axial.
 * 		θ^(bb)_(4i-2)=ArcCos[({x_(4i-3)-x_(4i-2)}∙{x_(4i-3)-x↓_(4i-2)})
 * 			/{‖x_(4i-3)-x_(4i-2)‖‖x_(4i-3)-x↓_(4i-2)‖}]
 * Por otro lado, los ángulos centrados en las posiciones 3, 7 y 11, se
 * encuentran formados por lo átomo anterior, el átomo actual y el átomo
 * anterior, otra vez, pero en la celda superior respecto a la dirección
 * axial.
 * 		θ^(bb)_(4i-1)=ArcCos[({x_(4i-1)-x_(4i-2)}∙{x_(4i-1)-x↑_(4i-2)})
 * 			/{‖x_(4i-1)-x_(4i-2)‖‖x_(4i-1)-x↑_(4i-2)‖}]
 * Y finalmente, los ángulos centrados en los átomos con posiciones 4, 8
 * y 12 se encuentran formados por el átomo siguiente, el átomo actual y
 * el átomo siguiente, ahora de la celda superior respecto a la
 * dirección axial.
 * 		θ^(bb)_4i=ArcCos[({x_4i-x_(4i+1)}∙{x_4i-x↑_(4i+1)})
 * 			/{‖x_4i-x_(4i+1)‖‖x_4i-x↑_(4i+1)‖}]
 *
 * 		Alto de celda. Sencillamente es la longitud del vector en la
 * dirección 3. Este parámetro no se obtiene por cada átomo, sino de
 * manera global.
 * 			z=‖a_3‖
 * 
 * 
 * 
 * ****Instalación y Uso.****
 * 
 * Para compilar el paquete basta con realizar la siguiente instrucción
 * en el directorio donde se encuentre el código fuente:
 * 
 * 		g++ -o cparam -Wall cparam.cpp -static-libgcc -static-libstdc++
 * 
 * lo que generará el archivo ejecutable cparam. En algunas ocasiones es
 * necesario cambiarle los privilegios a ejecutable mediante la
 * instrucción
 * 
 * 		chmod u+x cparam
 * 
 * Ya con ello, se podrá ejecutar por cualquier usuario. Siempre y
 * cuando se especifique la ruta completa.
 * 
 * Si no se desea especificar la ruta cada vez, puede agregarse a la
 * variable de entorno PATH, añadiendo el siguiente comando en el
 * archivo ".bashrc" del usuario
 * 
 * 		export PATH=$PATH:[ruta de instalación]/ 
 * 
 * Donde [ruta de instalación] es la ruta de donde se realizó la
 * instalación, sin poner corchetes y empleando el carácter “/” como
 * separador de directorios.
 * 
 * Por otro lado, la ejecución del paquete cparam es la siguiente:
 * 
 * 		cparam [nombre] [rutaCSV]
 * 
 * donde nombre indica el nombre del archivo a cargar. De no
 * especificarse se considerará que se desea cargar el archivo de salida
 * de VASP tipo CONTCAR en el directorio de trabajo. Se leen las bases y
 * se arroja toda la información de cada cálculo realizada, incluyendo
 * las coordenadas cartesianas, así como cada ángulo y longitud de
 * enlace individual. Finalmente, se calcularán las propiedades promedio
 * y se reportarán en pantalla. EL paquete igual permite el uso de un
 * segundo parámetro, una ruta para exportar las coordenadas en formato
 * "*.csv", Este parámetro será leído siempre y cuando se haya
 * especificado el parámetro nombre
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

using namespace std;

float longitud(float *v){
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

float distancia(float *v1, float *v2){
	float x = v1[0]-v2[0];
	x *= x;
	float y = v1[1]-v2[1];
	y *= y;
	float z = v1[2]-v2[2];
	z *= z;

	return sqrt(x+y+z);
}

float distanciacentral(float *v1, float *v2){
	float x = 0.5*(v1[0]+v2[0]);
	x *= x;
	float y = 0.5*(v1[1]+v2[1]);
	y *= y;
	float z = 0.5*(v1[2]+v2[2]);
	z *= z;
	return sqrt(x+y+z);
}

float radio(float *v){
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

float angulo(float *v1, float *v2, float *v3){
	float ax = v1[0] - v2[0], ay = v1[1] - v2[1], az = v1[2] - v2[2];
	float bx = v3[0] - v2[0], by = v3[1] - v2[1], bz = v3[2] - v2[2];

	return acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax + ay*ay + az*az)*(bx*bx + by*by + bz*bz)));
}

float angulosimetricoz(float *v1, float *v2){
	float ax = v1[0] - v2[0], ay = v1[1] - v2[1], az = v1[2] - v2[2];
	float axx = ax*ax, ayy = ay*ay, azz = az*az;

	return acos((axx+ayy-azz)/(axx + ayy + azz));
}

int main(int argc, char **args) {
	ifstream coordenadas;
	float factor;
	float bases[9];
	string comodin;
	
	float **atomos;
	float *coordth; //coordenada teta para cada átomo

	float **escalados;
	float **superiores;
	float **inferiores;
	
	string cadena; //cadena auxiliar

	int especies = 0,
	    num = 0, //número de átomos
	    numero;
	
	float n, ra = 0., rb = 0., thab= 0., thbb= 0., z, A = 0., D = 0.;

	cout << std::setprecision(17);

	if(argc >=2){
		cout <<"\n\nArchivo de entrada '" <<args[1] <<"'\n\n";
		string str(args[1]);
		coordenadas.open(str);
	}
	else{
		cout <<"\n\nNo se asignó un archivo de entrada, se empleará 'CONTCAR' para dicho fin\n\n";
		coordenadas.open("CONTCAR");
	}


	getline(coordenadas, cadena); //nombre del trabajo

	coordenadas //>>comodin >>comodin
		>>factor;

	cout <<"factor : " <<factor <<"\n";

	for(int i = 0; i<9; i++){
		coordenadas >>bases[i];
		bases[i] *= factor;
	}

	do{getline(coordenadas, cadena);} //nombre de especies
	while(cadena.find_first_not_of (' ') == cadena.npos && cadena.find_first_not_of ('\t') == cadena.npos);
	{
		stringstream linea(cadena);

		while(linea >> comodin){
			especies++;
		}
	}
	
	do{getline(coordenadas, cadena);} //atomos por especie
	while(cadena.find_first_not_of (' ') == cadena.npos && cadena.find_first_not_of ('\t') == cadena.npos);
	{
		stringstream linea(cadena);

		while(linea >> numero){
			num += numero;
		}
	}
	//coordenadas >>num; //átomos por especie
	coordenadas >>comodin;//donde se indica el tipo de coordenadas. Ahora se asume Direct
	
	n = num/4;
	
	atomos = new float*[num];
	escalados = new float*[num];
	inferiores = new float*[num];
	superiores = new float*[num];

	cout <<"Salidas escaladas\n";

	float centro[]={0., 0., 0.};
	
	for(int i = 0; i<num; i++){
		atomos[i] = new float[3];
		escalados[i] = new float[3];
		inferiores[i] = new float[3];
		superiores[i]= new float[3];

		coordenadas >>atomos[i][0] >>atomos[i][1] >>atomos[i][2];
		
		if(atomos[i][0] > 0.75)//experimental
			atomos[i][0] -=1.;
		
		if(atomos[i][1] > 0.75)//experimental
			atomos[i][1] -=1.;
		
		escalados[i][0] = bases[0]*atomos[i][0] + bases[3]*atomos[i][1] + bases[6]*atomos[i][2];
		escalados[i][1] = bases[1]*atomos[i][0] + bases[4]*atomos[i][1] + bases[7]*atomos[i][2];
		escalados[i][2] = bases[2]*atomos[i][0] + bases[5]*atomos[i][1] + bases[8]*atomos[i][2];
		
		centro[0] += escalados[i][0];
		centro[1] += escalados[i][1];
		centro[2] += escalados[i][2];

		cout <<"\t" <<escalados[i][0]<<", " <<escalados[i][1] <<", " <<escalados[i][2] <<"\n";
	}
	
	centro[0] /= num;
	centro[1] /= num;
	centro[2] /= num;

	coordth = new float[num];

	for(int i = 0; i < num; i++){
		escalados[i][0] -= centro[0];
		escalados[i][1] -= centro[1];
		escalados[i][2] -= centro[2];

		//angulo para posteriormente ordenar
		coordth[i] = atan2(escalados[i][1], escalados[i][0]);
	}
	
	/*for(int i = 0; i < num; i++){
		cout <<"\nth" <<i+1 << " = " <<coordth[i];
	}*/

	//ordenado de angulos
	for(int i = 0; i < num/2; i++){
		int maxi = num-1-i;
		int mini = i;
		float minth = coordth[mini];
		float maxth = coordth[maxi];

		if(minth > maxth){
			coordth[mini] = maxth;
			coordth[maxi] = minth;
			maxth = minth;
			minth = coordth[mini];

			float *tempcoord = escalados[mini];
			escalados[mini] = escalados[maxi];
			escalados[maxi] = tempcoord;
		}

		for(int j = i + 1; j< num -1 - i; j++){
			if(coordth[j] > maxth){
				coordth[maxi] = coordth[j];
				coordth[j] = maxth;
				maxth = coordth[maxi];

				float *tempcoord = escalados[j];
				escalados[j] = escalados[maxi];
				escalados[maxi] = tempcoord;
			}
			if(coordth[j] < minth){
				coordth[mini] = coordth[j];
				coordth[j] = minth;
				minth = coordth[mini];

				float *tempcoord = escalados[mini];
				escalados[mini] = escalados[j];
				escalados[j] = tempcoord;
			}
		}
	}

	//seleccion del ajuste de corrimiento de los indices
	int ajuste = 0;
	if(escalados[num-2][2] < 0. && escalados[num-1][2] < 0.){
		ajuste = num - 2;
	}
	if(escalados[num-1][2] < 0. && escalados[0][2] < 0.){
		ajuste = num - 1;
	}
	if(escalados[1][2] < 0. && escalados[2][2] < 0.){
		ajuste = 1;
	}

	//aplicacion del ajuste
	if(ajuste != 0){
		float **escaladostemp = new float*[num];
		float *coordthtemp = new float[num];

		for(int i = 0; i< num; i++){
			escaladostemp[i] = escalados[(i+ajuste) % num];
			coordthtemp[i] = coordth[(i+ajuste) % num];
		}

		delete[] escalados;
		delete[] coordth;
		escalados = escaladostemp;
		coordth = coordthtemp;
	}
	
	//impresion de las coordenadas ordenadas
	for(int i = 0; i < num; i++){
		cout <<"\t" <<escalados[i][0]<<", " <<escalados[i][1] <<", " <<escalados[i][2] <<"\n";
	}

	/*for(int i = 0; i < num; i++){
		cout <<"\nth" <<i+1 << " = " <<coordth[i];
	}*/

	//crea los espejos inferiores y superiores de los atomso en funcion de la direccion transversal de la celda unitaria.
	for(int i = 0; i < num; i++){
		inferiores[i][0] = escalados[i][0]-bases[6];
		inferiores[i][1] = escalados[i][1]-bases[7];
		inferiores[i][2] = escalados[i][2]-bases[8];

		superiores[i][0] = escalados[i][0]+bases[6];
		superiores[i][1] = escalados[i][1]+bases[7];
		superiores[i][2] = escalados[i][2]+bases[8];
	}

	//exportacion de las coordenadas cartesianas en formato *.csv
	if(argc >=3){
		//cout <<"\n\nArchivo de entrada '" <<args[1] <<"'\n\n";
		ofstream csv;
		string str(args[2]);
		csv.open(str);
		for(int i = 0; i < num; i++){
			csv <<escalados[i][0] <<"," <<escalados[i][1] <<"," <<escalados[i][2] <<"\n";
		}
		csv.close();
	}

	//reporta el centro
	cout <<"\n\nCentro : (" <<centro[0] <<", " <<centro[1] <<", " <<centro[2] <<")\n";

	coordenadas.close();
	
	//reporta las bases
	cout <<"Bases\n"
		<<"\t" <<bases[0] <<", " <<bases[1] <<", " <<bases[2] <<"\n"
		<<"\t" <<bases[3] <<", " <<bases[4] <<", " <<bases[5] <<"\n"
		<<"\t" <<bases[6] <<", " <<bases[7] <<", " <<bases[8] <<"\n";

	//reporta cada parametro medido
	cout <<"\nParámetros medidos\n";
	for(int i = 0; i < num; i += 2){
		float rai = distancia(escalados[i], escalados[i + 1]);
		float rbi = distancia(escalados[i + 1], escalados[(i + 2) % num]);
		float Ai = distanciacentral(escalados[i], escalados[i + 1]);
		//float thbbi1 = angulosimetricoz(escalados[(i+num-1) % num],escalados[i]);
		//float thbbi2 = angulosimetricoz(escalados[(i+2) % num],escalados[i+1]);
		//float thbbi1 = angulo(superiores[i],escalados[(i+num-1) % num],escalados[i]);
		//float thbbi2 = angulo(superiores[i+1],escalados[(i+2) % num],escalados[i+1]);
		float thbbi1 = angulo((i % 4 == 0?inferiores:superiores)[(i+num-1) % num],escalados[i],escalados[(i+num-1) % num]);
		float thbbi2 = angulo((i % 4 == 0?inferiores:superiores)[(i+2) % num],escalados[i+1],escalados[(i+2) % num]);
		//float thbbi1 = angulo(escalados[(i+num-1) % num],escalados[i],inferiores[(i+num-1) % num]);
		//float thbbi2 = angulo(escalados[(i+2) % num],escalados[i+1],inferiores[(i+2) % num]);
		float thabi1 = angulo(escalados[(i+num-2) % num],escalados[(i+num-1) % num],escalados[i]);
		float thabi2 = angulo(escalados[(i+3) % num],escalados[(i+2) % num],escalados[i+1]);
		float d1 = 2. * radio(escalados[i]);
		float d2 = 2. * radio(escalados[i + 1]);

	cout <<"\tr_a = " <<rai <<", r_b = " <<rbi <<", A = " <<Ai <<", th_ab = " <<thabi1 <<", th_ab = " <<thabi2 <<", th_bb = " <<thbbi1 <<", th_bb = " <<thbbi2 <<", D = " <<d1 <<", D = " <<d2 <<"\n";

		ra += rai;
		rb += rbi;
		A += Ai;
		thbb += thbbi1 + thbbi2;
		thab += thabi1 + thabi2;
		D += d1 +d2;
	}

	ra /= num / 2;
	rb /= num / 2;
	A /= num / 2;
	thbb /= num;
	thab /= num;
	D /= num;
	z = bases[8];

	//reporte de los promedios
	cout <<"\nParámetros promedios\n\tn = " <<n
		<<"\n\tra = " <<ra
		<<"\n\trb = " <<rb
		<<"\n\tA = " <<A
		<<"\n\tD = " <<D
		<<"\n\tU = " <<0.5*(longitud(bases)+longitud(bases+3))
		<<"\n\tth_ab = " <<thab
		<<"\n\tth_bb = " <<thbb
		<<"\n\tz = " <<z <<"\n";


	for(int i = 0; i<num; i++){
		delete[] atomos[i];
		delete[] escalados[i];
		delete[] inferiores[i];
		delete[] superiores[i];
	}
	delete[] coordth;
	delete[] atomos;
	delete[] escalados;
	delete[] inferiores;
	delete[] superiores;

	return 0;
}
