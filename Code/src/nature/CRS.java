package nature;

//Compressed Row Storage format for a sparse square matrix

public class CRS{

	//the values of the nonzero elements of the matrix
	float[] val;

	//the column indices of the elements in the val vector
	int[] col_idx;

	//locations in the val vector that start a row
	//the size is decided by the size of the matrix
	int[] row_ptr;

	//the number of rows of the matrix
	int mSize=0;


	//constructor that takes a sparse matrix and convert it to a CRS object
	CRS( float[][] matrix){

		int i;//row index
		int j;//column index

		int totalNonZeros;//the total number of nonzero in the matrix
		int index;

        //get the number of rows and columns
		mSize = matrix.length;
		//System.out.println("the size of the matrix is:"+mSize+"X"+mSize+"\n");

		//find the number of nonzero entries in the matrix
		totalNonZeros = 0;
		for(i=0; i<mSize; i++){//each row
			for(j=0; j<mSize; j++){//each column
				if(matrix[i][j] != 0)
					totalNonZeros++;
			}
		}

		//allocate memory for val and col_idx
		val = new float[totalNonZeros];
		col_idx = new int[totalNonZeros];

		//allocate memory for row_ptr
		row_ptr = new int[mSize+1];
		row_ptr[0] = 0;

		//store the matrix in CRS format
		index = 0;// point to the next position in val to store the value
		for(i=0; i<mSize; i++){

			for(j=0; j<mSize; j++){
				if(matrix[i][j] != 0){
					//add the value to val
					val[index] = matrix[i][j];
					//record the column index in col_idx
					col_idx[index] = j;
					index++;
				}
			}
			//update row_ptr
			row_ptr[i+1] = index;
		}
	}

	//print the matrix in CRS format
	public void printCRS(){
		int i;

		System.out.println("print vectors used in CRS");

		System.out.println("The vector val");
		for(i=0; i<val.length; i++){
			System.out.print( val[i] + ", ");
		}
		System.out.println();

		System.out.println("The vector col_idx");
		for(i=0; i<col_idx.length; i++){
			System.out.print( col_idx[i] + ", ");
		}
		System.out.println();

		System.out.println("The vector row_ptr");
		for(i=0; i<row_ptr.length; i++){
			System.out.print(row_ptr[i]+ ", ");
		}
		System.out.println();
	}

	//get the element in row i and column j in the matrix
	public float getElement(int i, int j){

		// Find where the row starts
		int rowStart = row_ptr[i];

		// Find where the next row starts
		int nextRowStarts = row_ptr[i+1];

		// Scan the column indices of entries in row i
		for(int k = rowStart; k < nextRowStarts; k++)
		{
			// if j is the column index of a non-zero, then return it
			if(col_idx[k] == j)
				return val[k];

			// if we have passed j, then entry (i, j) must be a zero
			if(col_idx[k] > j)
				return 0;
		}

		// if we have reached the end of the non-zeroes without
		// find j, it must be the index of a trailing 0
		return 0;

	}

	//print the matrix stored in CRS in a two dimensional array format
	public void printMatrix(){
		int i, j, k, zeroIndex;

		System.out.println("print the matrix in CRS format as a "+mSize+"x"+mSize+" matrix");

		//Scan the row_ptr array to find the beginning and end of row i
		for(i = 0; i < mSize; i++){

			//the column index of the element to print
			zeroIndex = 0;

			// Print row i: index k goes from the beginning to the end of row i
			for(k = row_ptr[i]; k < row_ptr[i+1]; k++){
				j = col_idx[k];

				//print entries of zero values that exist between consecutive
				// non-zeroes
				while(zeroIndex < j){
					System.out.print(0.0+", ");
					zeroIndex++;
				} // end of while-loop

				//print the nonzero value
				System.out.print(val[k]+", ");

				// Prepare zeroIndex for the next sequence of zeroes
				zeroIndex++;
			} // end of for-k loop

			//print the trailing zeroes in this row
			while(zeroIndex < mSize){
				System.out.print(0.0+", ");
				zeroIndex++;
			}

			// Start the next line
			System.out.println();
		} // end of for-i loop

		System.out.println();

	} // end of printMatrix


	/*takes a vector x and returns the product of the matrix stored in the CRS object with x.*/
	public float[] product(float[] x){
		//create vector to save product
		float[] product = new float[mSize];

		for (int i = 0; i < mSize; i++){
			for( int j = row_ptr[i]; j < row_ptr[i+1]; j++){
				product[i] += val[j] * x[col_idx[j]];

			}
		}

		return product;

	}//end of product

}




