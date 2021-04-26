

#pragma once

#include <utility>

namespace gadg{

    typedef std::pair<int,int> MatrixIndex;
    typedef std::tuple<Eigen::MatrixXi,Eigen::MatrixXi,Eigen::MatrixXi> SmithDecomposition;




    // Compute system A
    Eigen::MatrixXi systemA(const gadg::QuaternionScalar<int>& q){
        Eigen::MatrixXi A(4,6);
        int a = q.value.w;
        int b = q.value.x;
        int c = q.value.y;
        int d = q.value.z;

        A <<  b, c,  d, -b, -c, -d,
                a, -d, c, -a, -d,  c,
                d, a, -b,  d, -a, -b,
                -c, b,  a, -c,  b,  -a;

        return A;
    }


    MatrixIndex minAij(const Eigen::MatrixXi& A, const int s){
        int cols = A.cols();
        int rows = A.rows();

        // init min value with abs
        Eigen::MatrixXi Aabs = A.block(s,s,rows-s,cols-s).cwiseAbs();

//        std::cout << "Minaij::Aabs = \n"<<Aabs<<std::endl;
        MatrixIndex AijMin = std::make_pair(s,s);
        int valMin = Aabs.maxCoeff();

        for(int i=s ; i<rows ;++i){
            for(int j=s ; j<cols ; ++j){
                if(Aabs(i-s,j-s)!=0 && Aabs(i-s,j-s) <= valMin){
                    AijMin = std::make_pair(i,j);
                    valMin = Aabs(i-s,j-s);
                }
            }
        }



        return AijMin;
    }


    Eigen::MatrixXi Identity(const int n){
        return Eigen::MatrixXi::Identity(n,n);
    }

    void swapRows(Eigen::MatrixXi& A, const int i, const int j){
        A.row(i).swap(A.row(j));
    }

    void swapColumns(Eigen::MatrixXi& A, const int i, const int j){
        A.col(i).swap(A.col(j));
    }

    void addToRow(Eigen::MatrixXi& A,const int i, const int factor, const int j) {
        A.row(i) += factor*A.row(j);
    }

    void addToColumn(Eigen::MatrixXi& A,const int i, const int factor, const int j) {
        A.col(i) += factor*A.col(j);
    }

    void changeSignRow(Eigen::MatrixXi& A,const int i){
        A.row(i) *=-1;
    }

    void changeSignColumn(Eigen::MatrixXi& A,const int i){
        A.col(i) *=-1;
    }

    // simply check that the entry at place (s,s) in A is the only non-zero in its
    //s-th column under s and s-th row at the right of s
    bool isLone(const Eigen::MatrixXi& A,const int s) {
        for(int j=s+1 ;j<A.cols();++j){
            if(A(s,j)!=0)
                return false;
        }

        for(int i=s+1 ;i<A.rows();++i){
            if(A(i,s)!=0)
                return false;
        }
        return true;
    }

    // find an element which is not divisible by M(s,s)
    MatrixIndex getNextEntry(const Eigen::MatrixXi& A, const int s){
        for(int i=s+1;i<A.rows();++i){
            for(int j=s+1;j<A.cols();++j){
                if(A(i,j)%A(s,s) !=0)
                    return std::make_pair(i,j);
            }
        }
        return std::make_pair(-1,-1); // nothing found
    }



    // Smith Normal Form S of the M x N matrix A
    // see https://www.dlfer.xyz/post/2016-10-27-smith-normal-form/ for more details about the algo
    SmithDecomposition SmithDecompositionComputation(const Eigen::MatrixXi& A){
        // init
        Eigen::MatrixXi L= Identity(A.rows());
        Eigen::MatrixXi R= Identity(A.cols());

        Eigen::MatrixXi M = A;

        const int maxs = std::min(M.rows(),M.cols());

        for(int s=0 ; s<maxs ; ++s){
            std::cout <<"step "<<s+1 <<"/"<<maxs <<std::endl;
            std::cout << "M = \n"<<M <<"\n" << std::endl;

            while(!(isLone(M,s))){
                MatrixIndex Mminij =minAij(M,s);
                swapRows(M,s,Mminij.first);
                swapRows(L,s,Mminij.first);
                swapColumns(M,s,Mminij.second);
                swapColumns(R,s,Mminij.second);

                // add for lines
                for(int x=s+1 ; x<M.rows() ; ++x){
                    if(M(x,s) != 0){
                        int factor = M(x,s) / M(s,s);
                        addToRow(M,x,-factor,s);
                        addToRow(L,x,-factor,s);
                    }
                }

                // add for columns
                for(int x=s+1 ; x<M.cols() ; ++x){
                    if(M(s,x) != 0){
                        int factor = M(s,x) / M(s,s);
                        addToColumn(M,x,-factor,s);
                        addToColumn(R,x,-factor,s);
                    }
                }

                if(isLone(M,s)){
                    MatrixIndex indexEntry = getNextEntry(M,s);
                    if(indexEntry.first >= 0){ // /todo check
                        addToRow(M,s,1,indexEntry.first);
                        addToRow(L,s,1,indexEntry.first);
                    }else{
                        if(M(s,s)<0){
                            changeSignRow(M,s);
                            changeSignRow(L,s);
                        }
                    }
                }
            }
        }
//        std::cout << "U = \n"<<L<<std::endl;
//        std::cout << "V = \n"<<R<<std::endl;

        // check that L*M*V = A
        std::cout << "Uinv = "<<L.cast<float>().inverse() << std::endl;
        std::cout << "LMV = "<<(L*A*R) << std::endl;
        std::cout << "A = "<<A << std::endl;
    return std::make_tuple(M,L,R);

    }


    // get the rank obtained after a SNF decomposition
    int getRankSNF(const SmithDecomposition& snfOfA){
        Eigen::MatrixXi S = std::get<0>(snfOfA);
        int rank =0;
        for(int i=0;i<S.rows();++i){
            if(S(i,i)!=0){
                rank+=1;
            }
        }
        return rank;
    }

    // solve Ax=b through the Smith Normal form of A called SnfOfA
    Eigen::VectorXf solve(const SmithDecomposition& snfOfA, const Eigen::VectorXi& b, const int rank ){
        Eigen::MatrixXi S = std::get<0>(snfOfA);
        //compute the inverse of S = [D,0 ; 0,0]
        Eigen::MatrixXf Dinv = Eigen::MatrixXf::Zero(rank,rank);
        for(int i=0;i<rank;++i){
            Dinv(i,i) = 1.0f/(float)S(i,i);
        }

        // compute the product c=U*b
        Eigen::MatrixXi U = std::get<1>(snfOfA);
        std::cout << "U  =\n"<<U<<std::endl;
        Eigen::VectorXi c = U*b;

        // finally compute the product Dinv*c
        Eigen::VectorXf y = Dinv*c.cast<float>();
        std::cout << "y  =\n"<<y<<std::endl;

        Eigen::MatrixXi V = std::get<2>(snfOfA);

        Eigen::VectorXf y_withZeros(V.cols());
        std::cout << "y with zeros before init =\n"<<y_withZeros<<std::endl;
        y_withZeros << y,Eigen::VectorXf::Zero(std::abs(V.cols()-U.cols()));

        std::cout << "y with zeros =\n"<<y_withZeros<<std::endl;

        Eigen::VectorXf x = V.cast<float>()*y_withZeros;
        std::cout << "solution =\n"<<x<<std::endl;

        return x;
    }

    // solve Ax=b through the Smith Normal form of A called SnfOfA
    bool hasIntegerSolutions(const Eigen::MatrixXf& Dinv,const Eigen::MatrixXi& U, const Eigen::VectorXi& b, const int rank ){

        Eigen::VectorXi c = U*b;

        // finally compute the product Dinv*c
        Eigen::VectorXf y = Dinv*c.cast<float>();

        bool hasIntegerSolutions;
        for(int i=0 ; i<y.rows();++i){
            if(std::fmod(y(i), static_cast<float>(1.0)) != 0.0)
                return false;
        }
        return true;

    }



}

