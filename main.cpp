#include "Weibull.cpp"
#include <iostream>

using namespace std;

int main()
{
    Weibull<float> test;

    test.Wb_Name="Testname";
    test.Wb_Unity="Hrs";
    test.Wb_Gamma=0; // start life

    vector<float> temptbf;

    temptbf.push_back(1806);
	temptbf.push_back(2079);
    temptbf.push_back(77);
	temptbf.push_back(158);
	temptbf.push_back(244);
	temptbf.push_back(335);
	temptbf.push_back(437);
	temptbf.push_back(535);
	temptbf.push_back(646);
	temptbf.push_back(766);
	temptbf.push_back(897);
	temptbf.push_back(1040);
	temptbf.push_back(1574);
	temptbf.push_back(2414);
	temptbf.push_back(2846);
	temptbf.push_back(3454);
	temptbf.push_back(4494);
	temptbf.push_back(1198);
	temptbf.push_back(1374);
	temptbf.push_back(1374);

    test.generate(temptbf);

    for(auto i=0u;i<test.getFti().size();i++)
        cout << "TBF | " <<test.getTBF(i)<< " | Fi | "<<test.getFti(i)<<endl;

    cout << "\nMTBF: " << test.Wb_MTBF     <<" "<<test.Wb_Unity<<endl;
    cout << "Sigma: "   << test.Wb_Sigma <<" "<<test.Wb_Unity<<endl;
    cout << "Eta: "    << test.Wb_Eta      <<" "<<test.Wb_Unity<<endl;
    cout << "Beta: "   << test.Wb_Beta     <<" "<<endl;

    cout <<"for 50%: "  <<test.getTime(0.5)<<" "<<test.Wb_Unity<<endl;
    cout <<"for 95%: "  <<test.getTime(0.95)<<" "<<test.Wb_Unity<<endl;
    cout <<"for 99%: "  <<test.getTime(0.99)<<" "<<test.Wb_Unity<<endl;
    return 0;
}
