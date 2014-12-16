#include <iostream>
#include<vector>
#include<set>
#include<map>
#include<queue>
#include<stack>
#include<string>
#include<algorithm>
#include<functional>
#include<iomanip>
#include<cstdio>
#include<cmath>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<climits>
#include <fstream>
#include <utility>
#include <time.h>
#include<stdio.h>
 using namespace std;
 int main(){
 /*int ip,op,examples;
	ifstream infile1,infile2;
	infile1.open("test_data");
	//infile2.open("output5");
	ofstream outfile2,outfile1;
	outfile1.open("test_data");*/
	freopen("data_for_class1_dos","r",stdin);
	freopen("data_for_class1_dos_raw","w",stdout);
	//int a,b;
	//scanf("%d%d",&a,&b);
	//printf("%d %d\n",a,b);





	int n=42,m=102;
	int i,j;
	for(i=1;i<=m;i++) {
		for(j=1;j<=n;j++) {
			double d;
			if(j==2||j==3||j==4||j==42) {
                string s;
                cin >> s;
                if (j == 2){
                    if (s.compare("icmp")==0){
                        cout << "3 ";
                    }else if (s.compare("tcp")==0)cout << "1 ";
                    else if (s.compare("udp")==0)cout << "2 ";
                }else if (j == 3){
                    if (s.compare("ecr_i")==0){
                        cout << "1 ";
                    }else if (s.compare("private")==0)cout << "2 ";
                    else if (s.compare("http")==0)cout << "3 ";
                    else if (s.compare("ftp_data")==0)cout << "4 ";
                    else if (s.compare("smtp")==0)cout << "5 ";
                    else if (s.compare("eco_i")==0)cout << "6 ";
                    else if (s.compare("finger")==0)cout << "7 ";
                    else if (s.compare("Aol")==0)cout << "8 ";
                    else if (s.compare("bgp")==0)cout << "9 ";
                    else if (s.compare("auth")==0)cout << "10 ";
                    else if (s.compare("courier")==0)cout << "11 ";
                    else if (s.compare("csnet_ns")==0)cout << "12 ";
                    else if (s.compare("ctf")==0)cout << "13 ";
                    else if (s.compare("daytime")==0)cout << "14 ";
                    else if (s.compare("discard")==0)cout << "15 ";
                    else if (s.compare("domain")==0)cout << "16 ";
                    else if (s.compare("domain_u")==0)cout << "17 ";
                    else if (s.compare("echo")==0)cout << "18 ";
                    else if (s.compare("efs")==0)cout << "19 ";
                    else if (s.compare("aexec")==0)cout << "20 ";
                    else if (s.compare("ftp")==0)cout << "21 ";
                    else if (s.compare("gopher")==0)cout << "22 ";
                    else if (s.compare("harvest")==0)cout << "23 ";
                    else if (s.compare("hostnames")==0)cout << "24 ";
                    else if (s.compare("http_2784")==0)cout << "25 ";
                    else if (s.compare("http_443")==0)cout << "26 ";
                    else if (s.compare("http_8001")==0)cout << "27 ";
                    else if (s.compare("imap4")==0)cout << "28 ";
                    else if (s.compare("IRC")==0)cout << "29 ";
                    else if (s.compare("iso_tsap")==0)cout << "30 ";
                    else if (s.compare("klogin")==0)cout << "31 ";
                    else if (s.compare("kshell")==0)cout << "32 ";
                    else if (s.compare("ldap")==0)cout << "33 ";
                    else if (s.compare("link")==0)cout << "34 ";
                    else if (s.compare("login")==0)cout << "35 ";
                    else if (s.compare("mtp")==0)cout << "36 ";
                    else if (s.compare("name")==0)cout << "37 ";
                    else if (s.compare("netbios_dgm")==0)cout << "38 ";
                    else if (s.compare("netbios_ns")==0)cout << "39 ";
                    else if (s.compare("netbios_ssn")==0)cout << "40 ";
                    else if (s.compare("netstat")==0)cout << "41 ";
                    else if (s.compare("nnsp")==0)cout << "42 ";
                    else if (s.compare("nntp")==0)cout << "43 ";
                    else if (s.compare("ntp_u")==0)cout << "44 ";
                    else if (s.compare("other")==0)cout << "45 ";
                    else if (s.compare("pm_dump")==0)cout << "46 ";
                    else if (s.compare("pop_2")==0)cout << "47 ";
                    else if (s.compare("pop_3")==0)cout << "48 ";
                    else if (s.compare("printer")==0)cout << "49 ";
                    else if (s.compare("red_i")==0)cout << "50 ";
                    else if (s.compare("remote_job")==0)cout << "51 ";
                    else if (s.compare("rje")==0)cout << "52 ";
                    else if (s.compare("shell")==0)cout << "53 ";
                    else if (s.compare("sql_net")==0)cout << "54 ";
                    else if (s.compare("ssh")==0)cout << "55 ";
                    else if (s.compare("sunrpc")==0)cout << "56 ";
                    else if (s.compare("supdup")==0)cout << "57 ";
                    else if (s.compare("systat")==0)cout << "58 ";
                    else if (s.compare("telnet")==0)cout << "59 ";
                    else if (s.compare("tftp_u")==0)cout << "60 ";
                    else if (s.compare("tim_i")==0)cout << "61 ";
                    else if (s.compare("time")==0)cout << "62 ";
                    else if (s.compare("urh_i")==0)cout << "63 ";
                    else if (s.compare("urp_i")==0)cout << "64 ";
                    else if (s.compare("uucp")==0)cout << "65 ";
                    else if (s.compare("uucp_path")==0)cout << "66 ";
                    else if (s.compare("vmnet")==0)cout << "67 ";
                    else if (s.compare("whois")==0)cout << "68 ";
                    else if (s.compare("X11")==0)cout << "69 ";
                    else if (s.compare("Z39_50")==0)cout << "70 ";






                }else if (j == 4){
                    if (s.compare("SF")==0){
                        cout << "1 ";
                    }else if (s.compare("RSTR")==0)cout << "2 ";
                    else if (s.compare("REJ")==0)cout << "3 ";
                    else if (s.compare("S0")==0)cout << "4 ";
                    else if (s.compare("S1")==0)cout << "5 ";
                    else if (s.compare("S2")==0)cout << "6 ";
                    else if (s.compare("S3")==0)cout << "7 ";
                    else if (s.compare("SH")==0)cout << "8 ";
                    else if (s.compare("OTH")==0)cout << "9 ";
                    else if (s.compare("RSTO")==0)cout << "10 ";
                    else if (s.compare("RSTOS0")==0)cout << "11 ";





                }else if (j == 42){
                    if (s.compare("normal.")==0){
                        cout << "0 ";
                    }else cout << "1 ";
                }
				//printf("%s ",s);
			}else{
				scanf("%lf",&d);
			}
			switch(j) {
				case 1:printf("%lf  ",d);break;
				case 2:break;
				case 3:break;
				case 4:break;
				case 5:printf("%lf  ",d);
					break;
				case 6:printf("%lf  ",d);
					break;
				case 7:printf("%lf  ",d);
					break;
				case 8:printf("%lf  ",d);
					break;
				case 9:printf("%lf  ",d);
					break;
				case 10:printf("%lf  ",d);
					break;
				case 11:printf("%lf  ",d);
					break;
				case 12:printf("%lf  ",d);
					break;
				case 13:printf("%lf  ",d);
					break;
				case 14:printf("%lf  ",d);
					break;
				case 15:printf("%lf  ",d);
					break;
				case 16:printf("%lf  ",d);
					break;
				case 17:printf("%lf  ",d);
					break;
				case 18:printf("%lf  ",d);
					break;
				case 19:printf("%lf  ",d);
					break;
                case 20:printf("%lf  ",d);break;
				case 21:printf("%lf  ",d);break;
				case 22:printf("%lf  ",d);break;
				case 23:printf("%lf  ",d);
					break;
				case 24:printf("%lf  ",d);
					break;
				case 25:printf("%lf  ",d);break;
				case 26:printf("%lf  ",d);break;
				case 27:printf("%lf  ",d);break;
				case 28:printf("%lf  ",d);
					break;
				case 29:printf("%lf  ",d);
					break;
				case 30:printf("%lf  ",d);
					break;
				case 31:printf("%lf  ",d);break;
				case 32:printf("%lf  ",d);break;
				case 33:printf("%lf  ",d);break;
				case 34:printf("%lf  ",d);break;
				case 35:printf("%lf  ",d);break;
				case 36:printf("%lf  ",d);
					break;
				case 37:printf("%lf  ",d);break;
				case 38:printf("%lf  ",d);break;
				case 39:printf("%lf  ",d);
					break;
				case 40:printf("%lf  ",d);break;
				case 41:printf("%lf  ",d);break;
				case 42:break;
				case 43:break;
			}
		}
		printf("\n");
	}
	return 0;
 }



