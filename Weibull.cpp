#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#include <map>


template<class T> class Weibull
{
public:
    //constructor/destructor
    Weibull(void);
    ~Weibull(void);

    //setter
    void setTBF(std::vector<T> const &);
    void setTBF(T const & TBF);

    //getter
    T getProb(T const &)const;
    T getTime(T const &) const;
    T getFti(unsigned int const &) const;
    T getTBF(unsigned int const &) const;
    std::vector<T> getTBF(void) const;
    std::vector<T> getFti(void) const;
    std::string ExprReg(void) const;

    //over
    void clear(void);
    void generate(void); //if use setTBF
    void generate(std::vector<T> const &); // ifn't use setTBF

    //attribut
    T Wb_MTBF;
    T Wb_Sigma;
    T Wb_Beta;
    T Wb_Eta;
    T Wb_Gamma;
    std::string Wb_Unity;
    std::string Wb_Name;


private:
    //constante of linear equation y=Ax+B
    T Wb_A;
    T Wb_B;
    std::vector<T> Wb_TBF;
    std::vector<T> Wb_Fti;

    //calcul
    T ExponentialRegress(void);
};

//constructor/destructor
template<class T> Weibull<T>::Weibull(void)
{
    this->Wb_MTBF=0.0f;
    this->Wb_Sigma=0.0f;
    this->Wb_Beta=0.0f;
    this->Wb_Eta=0.0f;
    this->Wb_Gamma=0.0f;
    this->Wb_Unity="None";
    this->Wb_Name="None";

    this->Wb_A=0.0f;
    this->Wb_B=0.0f;
}

template<class T> Weibull<T>::~Weibull(void)
{
}

using std::vector;
using std::sort;

//Setter
template<class T>void Weibull<T>::setTBF(vector<T> const & dataTBF)
{
    this->Wb_TBF.clear();
    this->Wb_TBF=dataTBF;
}
template<class T>void Weibull<T>::setTBF(T const & TBF)
{
    this->Wb_TBF.push_back(TBF);
}
//Getter
template<class T> T Weibull<T>::getProb(T const & Time)const
{

    return exp(-pow((Time-this->Wb_Gamma)/this->Wb_Eta,this->Wb_Beta));
}

template<class T> T Weibull<T>::getTime(T const & Prob) const
{
    return this->Wb_Eta*pow(log(1/Prob),(1/this->Wb_Beta))+this->Wb_Gamma;
}

template<class T> T Weibull<T>::getFti(unsigned int const & Idx) const
{
    return this->Wb_Fti[Idx];
}

template<class T> T Weibull<T>::getTBF(unsigned int const & Idx) const
{
    return this->Wb_TBF[Idx];
}

template<class T> vector<T> Weibull<T>::getFti(void) const
{
    return this->Wb_Fti;
}

template<class T> vector<T> Weibull<T>::getTBF(void) const
{
    return this->Wb_TBF;
}

//Over
template<class T> void Weibull<T>::clear(void)
{
    this->Wb_MTBF=0;
    this->Wb_Sigma=0;
    this->Wb_Beta=0;
    this->Wb_Eta=0;
    this->Wb_Gamma=0;
    this->Wb_Unity="None";
    this->Wb_Name="None";

    this->Wb_A=0;
    this->Wb_B=0;

    this->Wb_TBF.clear();
    this->Wb_Fti.clear();
}

template<class T> void Weibull<T>::generate(void)
{
//calcule MTBF
    this->Wb_MTBF=0;
    for(auto it:this->Wb_TBF)
        this->Wb_MTBF+=it;
    this->Wb_MTBF=this->Wb_MTBF/this->Wb_TBF.size();

//Mise en orgdre croissant
    sort(this->Wb_TBF.begin(),this->Wb_TBF.end());

//calcule Fi
    this->Wb_Fti.clear();
    while(this->Wb_Fti.size()<this->Wb_TBF.size())
    {
        if(this->Wb_TBF.size()<=20)
            this->Wb_Fti.push_back(static_cast<T>(this->Wb_Fti.size()+1-0.3)/static_cast<T>(this->Wb_TBF.size()+0.4)); //fi=(i-0.3)/(n+0.4)

        else if(this->Wb_TBF.size()>20 && this->Wb_TBF.size()<=50)
            this->Wb_Fti.push_back((static_cast<T>(this->Wb_Fti.size()+1)/static_cast<T>(this->Wb_TBF.size()+1))); //fi=(i)/(n+1)

        else if(this->Wb_TBF.size()>50)
            this->Wb_Fti.push_back(static_cast<T>(this->Wb_Fti.size()+1)/static_cast<T>(this->Wb_TBF.size())); //fi=i/n
    }

//calcule ecart type
    for(auto i=0u;i<this->Wb_TBF.size();i++)
        this->Wb_Sigma+=(this->Wb_TBF[i]-this->Wb_MTBF)*(this->Wb_TBF[i]-this->Wb_MTBF);

    this->Wb_Sigma=sqrt(this->Wb_Sigma/static_cast<T>(this->Wb_TBF.size()));//E(x)=sqrt((1/n)(x*x+x*x+x*x,...)-moyenne*moyenne)

//regression lineaire

    this->ExponentialRegress();
//calcule du parametre d'echelle
    this->Wb_Eta=exp(log((1-exp(-1))/this->Wb_B)/this->Wb_A); //log((1-exp(-1))/b)/a=log(t)

//calcule du parametre de pente (derivée de at+b) etant dinne que par equivalance on a y=_b(t-ln(n)) ou y=ln(-ln(r(t)))
    this->Wb_Beta=-log(-log(1-(this->Wb_B*exp(this->Wb_A*(log(this->Wb_Eta)-1)))));
}

template<class T> void Weibull<T>::generate(vector<T> const & dataTBF)
{
    this->setTBF(dataTBF);
    this->generate();
}

//courbe tendance
template<class T> T Weibull<T>::ExponentialRegress(void)
{
    T xsomme{0}, ysomme{0}, xysomme{0}, xxsomme{0};

    for (auto i=0u;i<this->Wb_TBF.size();i++)
    {
        xsomme+=log(this->Wb_TBF[i]);
        ysomme+=log(this->Wb_Fti[i]);
        xysomme+=log(this->Wb_TBF[i])*log(this->Wb_Fti[i]);
        xxsomme+=log(this->Wb_TBF[i])*log(this->Wb_TBF[i]);
    }

    this->Wb_A = (static_cast<T>(this->Wb_TBF.size())*xysomme - xsomme*ysomme)/(static_cast<T>(this->Wb_TBF.size())*xxsomme - xsomme*xsomme);
    this->Wb_B = exp((ysomme - this->Wb_A*xsomme)/static_cast<T>(this->Wb_TBF.size()));

    return 0;//coefficient de dermination
}
