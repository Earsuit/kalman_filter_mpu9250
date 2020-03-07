/*
  Matrix.h - A free library for matrix calculation. Invalid dimension for calculation
  will lead to empty matrix, therefore, one should check whether the matrix is
  empty after calculation.

  This library is easy to use, I have overloaded operators, one can use +,*,/,-,
  =,==,!= as they want, but there is no definition for matrix division (like A/B).
  One should use the inverse function inv() instead.

  The transpose() and inv() functions are static functions, which should use
  Matrix<Type>::transpose() and Matrix<Type>::inv() format to call, Type is the
  type you used to create the matrix.

  The elements of the matrix are easy to access. For example, we have a 2 by 1 vector called A, one could use 
  A[1][0] to get the second element or A[1][0] to assign a value to it.

  The show(int decimal) function can be used to visulise the matrix, where decimal specifies the decimal numbers.

  Detailed description: https://playground.arduino.cc/Code/Matrix.

  Feedback and contribution is welcome!

  Version 1.3
  * Overloading operator [], so we can access the element by using A[][]. Because the internal data structure is a two 
    dimensional array, we have to use A[*][0] or A[0][*] to access the element in a column or row matrix.

  Version 1.2
  * Now it can perform mixed calculation with scalars. 
    But you need to ensure the dimension of the matrix is 1 for plus and minus.

  Version 1.1
  * Improved the processing speed, it becomes more efficient.

  ----------------
  Version 1.0

  Created by Yudi Ren, Jan 05, 2018.
  renyudicn@outlook.com
*/

#ifndef MATRIX_H
#define MATRIX_H

template <class Any>
class Matrix
{
    public:
        int _row;
        int _column;
        Any** _entity;
        Matrix();
        //ini stands for the initial value of all elements, r is row, c is column
        Matrix(int r,int c, Any ini=0);
        //copy constructor, one can directly use the "=" operator like A=B
        Matrix(const Matrix & m);
        //move constructor, improve the processing speed
        Matrix(Matrix && m);
        //for a 2d array, you need to cast it to (Any*) type, see example for
        //more information
        Matrix(int r,int c,Any* m);
        /*
            'I': create an identity matrix when r==c, otherwise it's empty
        */
        Matrix(int r,int c, char type);
        ~Matrix();
        Matrix & operator*=(Any a);
        Matrix & operator*=(const Matrix & A);
        Matrix & operator/=(Any a);
        Matrix & operator+=(const Matrix & A);
        Matrix & operator-=(const Matrix & A);
        Matrix operator*(const Matrix & A);
        Matrix operator*(Any a);
        template <class T>
        friend Matrix<T> operator*(T a, Matrix<T> & A);
        Matrix operator/(Any a);
        Matrix operator+(const Matrix & A);
        Any operator+(Any a);

        Any* operator[](int index);

        template <class T>
        friend T operator+(T a, const Matrix<T> & A);
        
        Matrix operator-(const Matrix & A);
        Any operator-(Any a);
        template <class T>
        friend T operator-(T a, const Matrix<T> & A);

        Matrix & operator=(const Matrix & A);
        //move assignment, improve the processing speed
        Matrix & operator=(Matrix && A);
        bool operator!=(const Matrix & A);
        bool operator==(const Matrix & A);
        //The inverse function uses the Gauss-Jordan Elimination,it won't modify
        //the original matrix and it returns the inverse of matrix A
        static Matrix inv(const Matrix & A);
        static Any inv(Any a);
        //The transpose function won't modify the original matrix and it returns
        //the transpose of matrix A
        static Matrix transpose(const Matrix & A);
        static Any transpose(Any a);
        //i: the row needs to be swaped, j: the position the row i goes to, this
        //function modifies the original matrix
        void swapRow(int i, int j);
        //decimal: the numbers of decimal place
        void show(int decimal = 0);
        //check if it's an empty matrix.
        bool notEmpty();
};

template <class Any>
Matrix<Any>::Matrix(){
    _row = 0;
    _column = 0;
    _entity = nullptr;
}

template <class Any>
Matrix<Any>::Matrix(int r,int c, Any ini){
    _row = r;
    _column = c;
    _entity = new Any*[_row];
    for(int i=0;i<_row;i++)
        _entity[i] = new Any[_column];

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++)
            _entity[i][j] = ini;
}

template <class Any>
Matrix<Any>::Matrix(const Matrix<Any> & m){
    _row = m._row;
    _column = m._column;

    _entity = new Any*[_row];
    for(int i=0;i<_row;i++)
        _entity[i] = new Any[_column];

    int bytes = _column*sizeof(Any);
    for(int i=0;i<_row;i++)
        memcpy(_entity[i],m._entity[i],bytes);
}

template <class Any>
Matrix<Any>::Matrix(Matrix<Any> && m){
    _row = m._row;
    _column = m._column;
    _entity = m._entity;
    m._entity = nullptr;
    m._row = 0;
    m._column = 0;
}

template <class Any>
Matrix<Any>::Matrix(int r,int c, Any* m){
    _row = r;
    _column = c;

    _entity = new Any*[_row];
    for(int i=0;i<r;i++)
        _entity[i] = new Any[_column];

    int bytes = _column*sizeof(Any);
    for(int i=0;i<_row;i++)
        memcpy(_entity[i],m+i*_column,bytes);

}

template <class Any>
Matrix<Any>::Matrix(int r,int c,char type){
    _row = r;
    _column = c;

    _entity = new Any*[_row];
    for(int i=0;i<r;i++)
        _entity[i] = new Any[_column];

    switch(type){
        case 'I':  //identity matrix
            if(r == c){
                for(int i=0;i<_row;i++)
                    for(int j=0;j<_column;j++)
                        _entity[i][j] = (i==j?1:0);
            }else{
                if(_entity != nullptr){
                    for(int i=0;i<_row;i++)
                        delete [] _entity[i];
                    delete [] _entity;
                }
            }
            break;
        default:
            for(int i=0;i<_row;i++)
                for(int j=0;j<_column;j++)
                    _entity[i][j] = 0;
    }
}

template <class Any>
Matrix<Any>::~Matrix(){
    if(_entity != nullptr){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];
        delete [] _entity;
    }
}

template <class Any>
Matrix<Any> Matrix<Any>::operator*(const Matrix<Any> & A){
    if(_column != A._row)
        Serial.println("Multiplication matrix dimension not match");

    Matrix<Any> tmp(_row,A._column);

    for(int i=0;i<_row;i++)
        for(int j=0;j<A._column;j++){
            for(int k=0;k<_column;k++)
                tmp._entity[i][j] += _entity[i][k]*A._entity[k][j];
        }

    return tmp;
}

template <class Any>
Matrix<Any> & Matrix<Any>::operator*=(const Matrix<Any> & A){
    //check if the dimension matches
    if(_column != A._row){
        if(_entity != nullptr){
            for(int i=0;i<_row;i++)
                delete [] _entity[i];
            delete [] _entity;
        }
        return *this;
    }

    Matrix<Any> tmp(*this);

    //the dimension of the result matrix does not match the calling matrix
    if(_column != A._column){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];

        _column = A._column;

        for(int i=0;i<_row;i++)
            _entity[i] = new Any[_column];

        for(int i=0;i<_row;i++)
            for(int j=0;j<_column;j++){
                _entity[i][j] = 0;
                for(int k=0;k<tmp._column;k++)
                    _entity[i][j] += tmp._entity[i][k] * A._entity[k][j];
            }
    }else{ //the dimension of the result matrix matches the calling matrix
        for(int i=0;i<_row;i++)
            for(int j=0;j<_column;j++){
                _entity[i][j] = 0;
                for(int k=0;k<_column;k++)
                    _entity[i][j] += tmp._entity[i][k]*A._entity[k][j];
            }
    }

    return *this;
}

template <class Any>
Matrix<Any> Matrix<Any>::operator*(Any a){
    Matrix<Any> tmp(_row,_column);

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            tmp._entity[i][j] = _entity[i][j]*a;
        }

    return tmp;
}

template <class Any>
Matrix<Any> & Matrix<Any>::operator*=(Any a){
    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            _entity[i][j] *= a;
        }

    return *this;
}

template <class Any>
Matrix<Any> operator*(Any a, Matrix<Any> & A){
    return A*a;
}

template <class Any>
Matrix<Any> Matrix<Any>::operator/(Any a){
    Matrix<Any> tmp(_row,_column);

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            tmp._entity[i][j] = _entity[i][j]/a;
        }

    return tmp;
}

template <class Any>
Matrix<Any> & Matrix<Any>::operator/=(Any a){
    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            _entity[i][j] /= a;
        }

    return *this;
}

template <class Any>
Matrix<Any> Matrix<Any>::operator+(const Matrix<Any> & A){
    if(_row != A._row || _column != A._column)
        Serial.println("Plus matrix dimension not match");

    Matrix<Any> tmp(_row,_column);

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            tmp._entity[i][j] = _entity[i][j]+A._entity[i][j];
        }

    return tmp;
}

template <class Any>
Any Matrix<Any>::operator+(Any a){
    if(_row != 1 || _column != 1)
        Serial.println("Plus scalar dimension not match");

    return _entity[0][0]+a;
}

template <class T>
T operator+(T a, const Matrix<T> & A){
    if(A._row != 1 || A._column != 1)
        Serial.println("Plus scalar dimension not match");
        
    return a+A._entity[0][0];
}

template <class Any>
Matrix<Any> & Matrix<Any>::operator+=(const Matrix<Any> & A){
    if(_row != A._row || _column != A._column){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];
        delete [] _entity;
        return *this;
    }

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            _entity[i][j] = _entity[i][j] + A._entity[i][j];
        }

    return *this;
}

template <class Any>
Matrix<Any> Matrix<Any>::operator-(const Matrix<Any> & A){
    if(_row != A._row || _column != A._column)
        Serial.println("Sub scalar dimension not match");

    Matrix<Any> tmp(_row,_column);

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            tmp._entity[i][j] = _entity[i][j]-A._entity[i][j];
        }

    return tmp;
}

template <class Any>
Any Matrix<Any>::operator-(Any a){
    if(_row != 1 || _column != 1)
        Serial.println("Sub scalar dimension not match");

    return _entity[0][0]-a;
}

template <class Any>
Any operator-(Any a, const Matrix<Any> & A){
    if(A._row != 1 || A._column != 1)
        Serial.println("Sub scalar dimension not match");
        
    return a-A._entity[0][0];
}

template <class Any>
Matrix<Any> & Matrix<Any>::operator-=(const Matrix<Any> & A){
    if(_row != A._row || _column != A._column){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];
        delete [] _entity;
        return *this;
    }

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++){
            _entity[i][j] = _entity[i][j] - A._entity[i][j];
        }

    return *this;
}

template <class Any>
Any* Matrix<Any>::operator[](int index){
    return _entity[index];
}

template <class Any>
Matrix<Any> Matrix<Any>::transpose(const Matrix<Any> & A){
    Matrix<Any> tmp(A._column,A._row);

    for(int i=0;i<A._row;i++)
        for(int j=0;j<A._column;j++)
            tmp._entity[j][i] = A._entity[i][j];

    return tmp;
}

template <class Any>
Matrix<Any>& Matrix<Any>::operator=(const Matrix<Any> & A){
    if(_entity != nullptr){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];
        delete [] _entity;
    }

    _row = A._row;
    _column = A._column;

    _entity = new Any*[_row];
    for(int i=0;i<_row;i++)
        _entity[i] = new Any[_column];

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++)
            _entity[i][j] = A._entity[i][j];

    return *this;
}

template <class Any>
Matrix<Any>& Matrix<Any>::operator=(Matrix<Any> && A){
    if(_entity != nullptr){
        for(int i=0;i<_row;i++)
            delete [] _entity[i];
        delete [] _entity;
    }
    _row = A._row;
    _column = A._column;
    _entity = A._entity;
    A._entity = nullptr;
    A._row = 0;
    A._column = 0;
    return *this;
}

template <class Any>
bool Matrix<Any>::operator!=(const Matrix<Any> & A){
    if(_row!=A._row || _column != A._column)
        return true;

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++)
            if(_entity[i][j] != A._entity[i][j])
                return true;

    return false;
}

template <class Any>
bool Matrix<Any>::operator==(const Matrix<Any> & A){
    if(_row != A._row || _column != A._column)
        return false;

    for(int i=0;i<_row;i++)
        for(int j=0;j<_column;j++)
            if(_entity[i][j] != A._entity[i][j])
                return false;

    return true;
}

template <class Any>
Matrix<Any> Matrix<Any>::inv(const Matrix & A){
    //has to be a square matrix
    if(A._row != A._column){
        Serial.println("Inv Has to be a square matrix");
    }

    Matrix<Any> tmp(A._row,A._column,'I');
    Matrix<Any> copy(A);

    //set the pivot to be the largest element in that column
    for(int i=0;i<A._row-1;i++){
        Any pivot = copy._entity[i][i];
        int row = i;
        for(int j = i+1;j<A._row;j++)
            if(copy._entity[j][i] > pivot){
                pivot = copy._entity[j][i];
                row = j;
            }

        if(pivot == 0){
            return Matrix<Any>();
        }else{
            copy.swapRow(i,row);
            tmp.swapRow(i,row);
        }
    }

    for(int i=0;i<A._row;i++)
        for(int j=0;j<A._row;j++)
            if(i!=j){
                Any backup = copy._entity[j][i];
                for(int k=0;k<A._column;k++){
                    copy._entity[j][k] -= backup*copy._entity[i][k]/copy._entity[i][i];
                    tmp._entity[j][k] -= backup*tmp._entity[i][k]/copy._entity[i][i];
                }
            }

    for(int i=0;i<A._row;i++)
        for(int j=0;j<A._column;j++){
            tmp._entity[i][j] /= copy._entity[i][i];
        }

    return tmp;
}

template <class Any>
Any Matrix<Any>::inv(Any a){
    return 1.0/a;
}

template <class Any>
void Matrix<Any>::swapRow(int i,int j){
    if(i!=j && i<_row && j<_row && i>=0 && j>=0){
        Any tmp;
        for(int k=0;k<_column;k++){
            tmp = _entity[j][k];
            _entity[j][k] = _entity[i][k];
            _entity[i][k] = tmp;
        }
    }
}

template <class Any>
void Matrix<Any>::show(int decimal)
{
    if(decimal<=0){
        for(int i=0; i<_row; i++)
    	{
    	    Serial.print("[");
    	    for(int j=0; j<_column; j++)
    		{
    			Serial.print("  ");
                Serial.print(_entity[i][j]);
                Serial.print(" ");
    		}
    		Serial.println("]");
    	}
        Serial.println();
    }else{
        for(int i=0; i<_row; i++)
    	{
    	    Serial.print("[");
    	    for(int j=0; j<_column; j++)
    		{
    			Serial.print("  ");
                Serial.print(_entity[i][j],decimal);
                Serial.print(" ");
    		}
    		Serial.println("]");
    	}
        Serial.println();
    }
}

template <class Any>
bool Matrix<Any>::notEmpty(){
    if(_entity == nullptr)
        return false;
    return true;
}
#endif
