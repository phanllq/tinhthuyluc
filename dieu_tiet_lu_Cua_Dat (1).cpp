#include <iostream>
#include <cmath>

using namespace std;
double Z[100],V[100],T[200],Qden[200],Qden0[200];
long n,m;
double epsilon=1,m_tran=0.35,Zng_tran=97,B=55,MNTL=110;
double qo,Qxa[200],V1[200],V2[200],Z1[200],Htr[200];
FILE *fptr;
double noi_suyZV(double ZZ)
{
	double VV=1789.4;
	long i=0;
	while ((Z[i]<ZZ)&&(i<=32))
		i++;
	VV=V[i-1]+(ZZ-Z[i-1])/(Z[i]-Z[i-1])*(V[i]-V[i-1]);
	return VV;
}
double noi_suyVZ(double VV)
{
	double ZZ=130;
	long i=0;
	while ((V[i]<VV)&&(i<=32))
		i++;
	ZZ=Z[i-1]+(VV-V[i-1])/(V[i]-V[i-1])*(Z[i]-Z[i-1]);
	return ZZ;
}
void nhap()
{
	FILE*f = fopen("lu_Cua_Dat.txt","rb");fptr = fopen("kq_0.5.out","w");
	fscanf(f,"%ld",&n);
	for(long i=0;i<n;i++)
	{
		fscanf(f,"%lf",&Z[i]);
		fscanf(f,"%lf",&V[i]);
	}
	for(long i=0;i<n;i++)
	{
		cout<<Z[i]<<" "<<V[i]<<endl;
	}
	fscanf(f,"%ld",&m);
	for(long i=0;i<m;i++)
	{
		fscanf(f,"%lf",&T[i]);
		fscanf(f,"%lf",&Qden0[i]);//m_tran=0.41;
		fscanf(f,"%lf",&Qden0[i]);//m_tran=0.32;
		fscanf(f,"%lf",&Qden[i]);m_tran=0.28;
		fscanf(f,"%lf",&Qden0[i]);//m_tran=0.24;
	}
	for(long i=0;i<m;i++)
	{
		cout<<T[i]<<" "<<Qden[i]<<endl;
	}

}

double tinhmax(double Bi,double MNTLi)
{
	B=Bi;
	MNTL=MNTLi;
	qo=epsilon*m_tran*B*pow((2*9.81),0.5)*pow((MNTL-Zng_tran),1.5);
	//cout<<"Luu luong xa qo "<<qo<<endl;
	bool dau_thoi_doan=true;
	for(long i=0;i<m;i++)
	{
		if (Qden[i]<=qo)
		{
			if (dau_thoi_doan)
			{
				Qxa[i]=Qden[i];
				Z1[i]=MNTL;
				V1[i]=noi_suyZV(Z1[i]);
				Htr[i]=MNTL-Zng_tran;
				V2[i]=V1[i];				
			}
			else
			{
				Qxa[i]=Qxa[i-1];
				V1[i]=V1[i-1]+((Qden[i-1]+Qden[i])/2-(Qxa[i-1]+Qxa[i])/2)*3600/pow(10,6);
				Htr[i]=noi_suyVZ(V1[i])-Zng_tran;
				double Qxa_moi=epsilon*m_tran*B*pow((2*9.81),0.5)*pow(Htr[i],1.5);
				int dem=0;
				while ((abs(Qxa[i]-Qxa_moi)>0.01)||(dem>100))
				{
					dem++;
					Qxa[i]=(Qxa[i]+Qxa_moi)/2;
					V1[i]=V1[i-1]+((Qden[i-1]+Qden[i])/2-(Qxa[i-1]+Qxa[i])/2)*3600/pow(10,6);
					Htr[i]=noi_suyVZ(V1[i])-Zng_tran;
					Qxa_moi=epsilon*m_tran*B*pow((2*9.81),0.5)*pow(Htr[i],1.5);
				}
				Z1[i]=Zng_tran+Htr[i];
				V2[i]=noi_suyZV(Z1[i]);
				if (Qxa[i]<=qo)
				{
					Qxa[i]=Qden[i];
					V1[i]=V1[i-1]+((Qden[i-1]+Qden[i])/2-(Qxa[i-1]+Qxa[i])/2)*3600/pow(10,6);
					Z1[i]=noi_suyVZ(V1[i]);
					Htr[i]=Z1[i]-Zng_tran;
				}				
			}
		}
		else
		{
			dau_thoi_doan=false;
			Qxa[i]=Qxa[i-1];
			V1[i]=V1[i-1]+((Qden[i-1]+Qden[i])/2-(Qxa[i-1]+Qxa[i])/2)*3600/pow(10,6);
			Htr[i]=noi_suyVZ(V1[i])-Zng_tran;
			double Qxa_moi=epsilon*m_tran*B*pow((2*9.81),0.5)*pow(Htr[i],1.5);
			int dem=0;
			while ((abs(Qxa[i]-Qxa_moi)>0.01)||(dem>100))
			{
				dem++;
				Qxa[i]=(Qxa[i]+Qxa_moi)/2;
				V1[i]=V1[i-1]+((Qden[i-1]+Qden[i])/2-(Qxa[i-1]+Qxa[i])/2)*3600/pow(10,6);
				Htr[i]=noi_suyVZ(V1[i])-Zng_tran;
				Qxa_moi=epsilon*m_tran*B*pow((2*9.81),0.5)*pow(Htr[i],1.5);
			}
			Z1[i]=Zng_tran+Htr[i];
			V2[i]=noi_suyZV(Z1[i]);
		}
		fprintf(fptr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",T[i],Qden[i],Qxa[i],Z1[i],V1[i],Htr[i],V2[i]);
		//cout<<T[i]<<"\t"<<Qden[i]<<"\t"<<Qxa[i]<<"\t"<<Z1[i]<<"\t"<<V1[i]<<"\t"<<Htr[i]<<"\t"<<V2[i]<<endl;
	}
	double max=Z1[0];
	for(long i=0;i<m;i++)
		if (max<Z1[i])
			max=Z1[i];
	return max;
}

int main()
{
	int NguoiDung=0;						// Menu lua chon nguoi dung
	cout<<"\n Hay lua chon : ";
	cout<<"\n 0 De thoat";
	cout<<"\n 1 Nguoi quan tri";
	cout<<"\n 2 Nguoi dung\n";
	cin>>NguoiDung;
	switch (NguoiDung) {
		case 1: {
			cout<<"\n Nguoi quan tri co the lua chon 5 va 6";
			break;
		}
		case 2: {
			cout<<"\n Nguoi dung co the lua chon tu 1 den 4";
			break;
		}
		default: 
			break;
	}
	int luaChon=0;
	do {
		if (NguoiDung>0) {
			cout<<"\n Hay lua chon : ";
			cout<<"\n 0 De thoat";
			if (NguoiDung==2) {
				cout<<"\n 1 Tinh kha nang xa - Q xa max";
				cout<<"\n 2 Xuat ket qua kha nang xa - Q xa max";
				cout<<"\n 3 Tinh dieu tiet - Q xa ~ t";
				cout<<"\n 4 Xuat ket qua - Q xa ~ t\n";
			}
			else if (NguoiDung==1) {
				cout<<"\n 5 Thong bao danh cho tan suat thiet ke va kiem tra";
				cout<<"\n 6 Cap nhat du lieu\n";
			}
			cin>>luaChon;
		}
		if (luaChon==0) break;
		else
			{
				double Nguong_moi,So_Khoang_moi,Pin_moi,Bien_moi,B_moi,epsilon_moi,qo_max;
				switch (luaChon) {
					case 1: {
						cout<<"\n + Nhap moi cao trinh nguong : ";
						cin>>Nguong_moi;
						cout<<" + Nhap moi so luong khoang tran : ";
						cin>>So_Khoang_moi;
						cout<<" + Nhap moi chieu day tru pin : ";
						cin>>Pin_moi;
						cout<<" + Nhap moi chieu day tru bien : ";
						cin>>Bien_moi;
						cout<<" + Nhap chieu rong mot khoang tran : ";
						cin>>B_moi;
						epsilon_moi=So_Khoang_moi*B_moi/(So_Khoang_moi*B_moi+(So_Khoang_moi-1)*Pin_moi);
						qo_max=epsilon_moi*m_tran*So_Khoang_moi*B_moi*pow((2*9.81),0.5)*pow((MNTL-Zng_tran),1.5);
						break;
					}
					case 2: {				
						cout<<"\nHe so co hep ben epsilon moi tinh duoc : "<<epsilon_moi<<endl;
						cout<<"Luu luong xa ung voi muc nuoc truoc lu tinh duoc : "<<qo_max<<endl;
						break;	
					}
					case 3: {				
						break;	
					}
					case 4: {				
						double B_tran[10],Z_truoc_lu[10],epsilon_theoB[10];
						B_tran[0]=55;
						B_tran[1]=44;
						B_tran[2]=33;
						B_tran[3]=30;
						B_tran[4]=20;
						B_tran[5]=10;
	
						epsilon_theoB[0]=0.846153846;
						epsilon_theoB[1]=0.854368932;
						epsilon_theoB[2]=0.868421053;
						epsilon_theoB[3]=0.882352941;
						epsilon_theoB[4]=0.909090909;
						epsilon_theoB[5]=1;

						Z_truoc_lu[0]=110;
						Z_truoc_lu[1]=22.1; // truong hop tinh nhieu Z
						Z_truoc_lu[2]=23.3;
						Z_truoc_lu[3]=24;
						Z_truoc_lu[4]=24.4;
						nhap();
	
						for(int Zi=0;Zi<1;Zi++)
						{
							for(int Bi=0;Bi<3;Bi++)
							{
								epsilon=epsilon_theoB[Bi];
								fprintf(fptr,"Truong hop B = %lf\t Muc nuoc ban dau = %lf\n",B_tran[Bi],Z_truoc_lu[Zi]);
								cout<<"Cac muc nuoc thuong luu lon nhat : "<<endl;
								cout<<tinhmax(B_tran[Bi],Z_truoc_lu[Zi])<<" ";			
							}
							cout<<endl;
						}
						fclose(fptr);
						break;	
					}
					case 5: {				
						break;	
					}
					case 6: {				
						break;	
					}
					default: 
						break;
			
				}
			}
		}
	while (luaChon>0);
	return 0;
}



