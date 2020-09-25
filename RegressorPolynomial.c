/*
	This file is part of RegressorPolynomial.

    RegressorPolynomial is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RegressorPolynomial is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RegressorPolynomial.  If not, see <https://www.gnu.org/licenses/>
*/

/*
	This program was developed by Eduardo Vasconcelos; associated professor on 
	Federal Institute of Science, Education and Technology of Pernambuco.
	
	If you have any question about this program or the Combinatorial Regression, 
	please, mail me at eduardo.vasconcelos@recife.ifpe.edu.br
*/

#include "csvReader.h"
#include "polynomialregression.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
* This program receives arguments, and the first one is obrigatory.
*  
* This program generates two files, the first is a .txt with the RSquared value
* and the regression formula in to be used in a spreadsheet tool,
* and the second is a .csv file with the residues obtined from the model.
*/
int main(int argc, char* argv[]){
	
	if(argc < 2){
		
		printf("Required arguments are missing! \n");
		printf("The correct form to start Regressor is: Regressor csvpath degree buffer_size \n");
		printf("You have to pass at last the csvpath argument \n");
		printf("if you don't define the degree argument, Regressor will consider the value 1 \n");
		printf("if you don't define the buffer_size, Regressor will consider the value 20'\n");
		
		exit(0);
	}
	
	int degree = 1;
	int buffer = 20;
	
	if(argc > 2){
	
		sscanf(argv[2], "%d", &degree);
	
		if(argc > 3)
			sscanf(argv[3], "%d", &buffer);
	}
	
	char* path = argv[1];
	
	double* X;
	double* Y;
	
	int rows;
	int variables;
	
	/*
	* Retrieving data from CSV file. X, y, rows and columns are passed as reference
	* parameters.
	*/
	getDataFromCSV(path, &X, &Y, &rows, &variables, buffer);
	
	/*
	* retrieving the beta coeficients
	*/
	double* B = performRegression(Y, X, degree, rows, variables);
    
    while (degree * variables > rows)
    {
        degree--;
    }
    
    int columns = degree * variables + 1;
    
    /*for(int i = 0; i < columns; i++){
    	printf("%lf \n",*(B + i));
	}*/

	/*
	* The following procedure is used to compute predicted array.
	* This array is used to compute the RSquared value and to genetate the 
	* residues.csv file.
	*/
    double* predicted = (double*) malloc(rows * sizeof(double));

	/*
	* This variable is used to calculate the total variation of array Y.
	*/
    double averageY = 0.0; 

    for (int i = 0; i < rows; i++)
    {
        averageY += *(Y + i);

		/*
		* Adding the beta_0 on polynomial
		*/
        *(predicted + i) = *B;

		/*
		* This code calculate each term of the polynomial.
		* The term structure is: beta_i x prod{X_i,j}
		*/
        for(int j = 1; j < columns; j++){

            *(predicted + i) += *(B + j) * pow(*(X + i*variables + ((j - 1) % variables)),((j - 1)/variables)+1);
            
        }
    }

    averageY /= rows;

    double sstot = 0.0;
    double ssres = 0.0;

    for (int i = 0; i < rows; i++)
    {
        sstot += pow((*(Y + i)) - averageY,2);
        ssres += pow((*(Y + i)) - *(predicted + i),2);
    }
    
    double RSquared = 1.0 - ssres/sstot;
    
    /*
    * Creating the spreadsheet formated formula and inserting it into a txt file.
    */
    printf("RSquared value is: %lf \n", RSquared);
    
    FILE * file;

    file = fopen("spreadsheet_equation.txt","w");

    fprintf(file, "this equation has a RSquared value of %lf \n", RSquared);

    fprintf(file, "put the following equation on cel A%i \n", variables+1);


    fprintf(file,"= %E + ",*B);

    for(int i = 1; i < columns; i++){
        
        fprintf(file,"(%E)*",*(B + i));

            int exp = (i + variables - 1)/variables; 

            if(exp == 1)
                fprintf(file, "(A%d)", ((i-1) % variables)+1);
            else
                fprintf(file, "(A%d^%d)",  ((i-1) % variables)+1, exp);

        if(i < columns - 1)
            fprintf(file, " + ");
        
    }

    fclose(file);
    
    /*
    * Creating the residues file.
    */
    file = fopen("residues.csv","w");
    
    fprintf(file, "Y,predicted,residuo\n");
    
    for(int i = 0; i < rows; i++){ 
    	fprintf(file, "%lf,%lf,%lf\n", *(Y + i), *(predicted + i), *(Y + i) - *(predicted + i));
	}
	
	fclose(file);
	
}
