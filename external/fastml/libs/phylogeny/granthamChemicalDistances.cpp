// $Id: granthamChemicalDistances.cpp 962 2006-11-07 15:13:34Z privmane $

#include "granthamChemicalDistances.h"
#include <cmath>

granthamChemicalDistances::granthamChemicalDistances() {
	for (int i=0; i<20;++i)		GranChemDist[i][i]=0;
	GranChemDist[0][1]=112;		GranChemDist[0][2]=111;		GranChemDist[0][3]=126;		GranChemDist[0][4]=195;		GranChemDist[0][5]=91;		GranChemDist[0][6]=107;
	GranChemDist[0][7]=60;		GranChemDist[0][8]=86;		GranChemDist[0][9]=94;		GranChemDist[0][10]=96;		GranChemDist[0][11]=106;	GranChemDist[0][12]=84;
	GranChemDist[0][13]=113;	GranChemDist[0][14]=27;		GranChemDist[0][15]=99;		GranChemDist[0][16]=58;		GranChemDist[0][17]=148;	GranChemDist[0][18]=112;
	GranChemDist[0][19]=64;
	
	GranChemDist[1][2]=86;		GranChemDist[1][3]=96;		GranChemDist[1][4]=180;		GranChemDist[1][5]=43;		GranChemDist[1][6]=54;		GranChemDist[1][7]=125;
	GranChemDist[1][8]=29;		GranChemDist[1][9]=97;		GranChemDist[1][10]=102;	GranChemDist[1][11]=26;		GranChemDist[1][12]=91;		GranChemDist[1][13]=97;
	GranChemDist[1][14]=103;	GranChemDist[1][15]=110;	GranChemDist[1][16]=71;		GranChemDist[1][17]=101;	GranChemDist[1][18]=77;		GranChemDist[1][19]=96;

	GranChemDist[2][3]=23;		GranChemDist[2][4]=139;		GranChemDist[2][5]=46;		GranChemDist[2][6]=42;		GranChemDist[2][7]=80;		GranChemDist[2][8]=68;
	GranChemDist[2][9]=149;		GranChemDist[2][10]=153;	GranChemDist[2][11]=94;		GranChemDist[2][12]=142;	GranChemDist[2][13]=158;	GranChemDist[2][14]=91;
	GranChemDist[2][15]=46;		GranChemDist[2][16]=65;		GranChemDist[2][17]=174;	GranChemDist[2][18]=143;	GranChemDist[2][19]=133;

	GranChemDist[3][4]=154;		GranChemDist[3][5]=61;		GranChemDist[3][6]=45;		GranChemDist[3][7]=94;		GranChemDist[3][8]=81;
	GranChemDist[3][9]=168;		GranChemDist[3][10]=172;	GranChemDist[3][11]=101;	GranChemDist[3][12]=160;	GranChemDist[3][13]=177;	GranChemDist[3][14]=108;
	GranChemDist[3][15]=65;		GranChemDist[3][16]=85;		GranChemDist[3][17]=181;	GranChemDist[3][18]=160;	GranChemDist[3][19]=152;

	GranChemDist[4][5]=154;		GranChemDist[4][6]=170;		GranChemDist[4][7]=159;		GranChemDist[4][8]=174;
	GranChemDist[4][9]=198;		GranChemDist[4][10]=198;	GranChemDist[4][11]=202;	GranChemDist[4][12]=196;	GranChemDist[4][13]=205;	GranChemDist[4][14]=169;
	GranChemDist[4][15]=112;	GranChemDist[4][16]=149;	GranChemDist[4][17]=215;	GranChemDist[4][18]=194;	GranChemDist[4][19]=192;

	GranChemDist[5][6]=29;		GranChemDist[5][7]=87;		GranChemDist[5][8]=24;
	GranChemDist[5][9]=109;		GranChemDist[5][10]=113;	GranChemDist[5][11]=53;		GranChemDist[5][12]=101;	GranChemDist[5][13]=116;	GranChemDist[5][14]=76;
	GranChemDist[5][15]=68;		GranChemDist[5][16]=42;		GranChemDist[5][17]=130;	GranChemDist[5][18]=99;		GranChemDist[5][19]=96;

	GranChemDist[6][7]=98;		GranChemDist[6][8]=40;
	GranChemDist[6][9]=134;		GranChemDist[6][10]=138;	GranChemDist[6][11]=56;		GranChemDist[6][12]=126;	GranChemDist[6][13]=140;	GranChemDist[6][14]=93;
	GranChemDist[6][15]=80;		GranChemDist[6][16]=65;		GranChemDist[6][17]=152;	GranChemDist[6][18]=122;	GranChemDist[6][19]=121;

	GranChemDist[7][8]=89;
	GranChemDist[7][9]=135;		GranChemDist[7][10]=138;	GranChemDist[7][11]=127;	GranChemDist[7][12]=127;	GranChemDist[7][13]=153;	GranChemDist[7][14]=42;
	GranChemDist[7][15]=56;		GranChemDist[7][16]=59;		GranChemDist[7][17]=184;	GranChemDist[7][18]=147;	GranChemDist[7][19]=109;

	GranChemDist[8][9]=94;		GranChemDist[8][10]=99;		GranChemDist[8][11]=32;		GranChemDist[8][12]=87;		GranChemDist[8][13]=100;	GranChemDist[8][14]=77;
	GranChemDist[8][15]=89;		GranChemDist[8][16]=47;		GranChemDist[8][17]=115;	GranChemDist[8][18]=83;		GranChemDist[8][19]=84;

	GranChemDist[9][10]=5;		GranChemDist[9][11]=102;	GranChemDist[9][12]=10;		GranChemDist[9][13]=21;		GranChemDist[9][14]=95;
	GranChemDist[9][15]=142;	GranChemDist[9][16]=89;		GranChemDist[9][17]=61;		GranChemDist[9][18]=33;		GranChemDist[9][19]=29;

	GranChemDist[10][11]=107;	GranChemDist[10][12]=15;	GranChemDist[10][13]=22;	GranChemDist[10][14]=98;
	GranChemDist[10][15]=145;	GranChemDist[10][16]=92;	GranChemDist[10][17]=61;	GranChemDist[10][18]=36;	GranChemDist[10][19]=32;

	GranChemDist[11][12]=95;	GranChemDist[11][13]=102;	GranChemDist[11][14]=103;
	GranChemDist[11][15]=121;	GranChemDist[11][16]=78;	GranChemDist[11][17]=110;	GranChemDist[11][18]=85;	GranChemDist[11][19]=97;

	GranChemDist[12][13]=28;	GranChemDist[12][14]=87;
	GranChemDist[12][15]=135;	GranChemDist[12][16]=81;	GranChemDist[12][17]=67;	GranChemDist[12][18]=36;	GranChemDist[12][19]=21;

	GranChemDist[13][14]=114;
	GranChemDist[13][15]=155;	GranChemDist[13][16]=103;	GranChemDist[13][17]=40;	GranChemDist[13][18]=22;	GranChemDist[13][19]=50;

	GranChemDist[14][15]=74;	GranChemDist[14][16]=38;	GranChemDist[14][17]=147;	GranChemDist[14][18]=110;	GranChemDist[14][19]=68;

	GranChemDist[15][16]=58;	GranChemDist[15][17]=177;	GranChemDist[15][18]=144;	GranChemDist[15][19]=124;

	GranChemDist[16][17]=128;	GranChemDist[16][18]=92;	GranChemDist[16][19]=69;

	GranChemDist[17][18]=37;	GranChemDist[17][19]=88;

	GranChemDist[18][19]=55;


	GranPolarityTable[0]=8.1  ; //A
    GranPolarityTable[1]=10.5 ; //R
    GranPolarityTable[2]=11.6 ; //N
    GranPolarityTable[3]=13.0 ; //D
    GranPolarityTable[4]=5.5  ; //C
    GranPolarityTable[5]=10.5 ; //Q
    GranPolarityTable[6]=12.3 ; //E
    GranPolarityTable[7]=9.0  ; //G
    GranPolarityTable[8]=10.4 ; //H
    GranPolarityTable[9]=5.2  ; //I
    GranPolarityTable[10]=4.9 ; //L
    GranPolarityTable[11]=11.3; //K
    GranPolarityTable[12]=5.7 ; //M
    GranPolarityTable[13]=5.2 ; //F
    GranPolarityTable[14]=8.0 ; //P
    GranPolarityTable[15]=9.2 ; //S
    GranPolarityTable[16]=8.6 ; //T
    GranPolarityTable[17]=5.4 ; //W
    GranPolarityTable[18]=6.2 ; //Y
    GranPolarityTable[19]=5.9 ; //V

/*
	GranVolumeTable[0]=8.1  ; //A
    GranVolumeTable[1]=10.5 ; //R
    GranVolumeTable[2]=11.6 ; //N
    GranVolumeTable[3]=13.0 ; //D
    GranVolumeTable[4]=5.5  ; //C
    GranVolumeTable[5]=10.5 ; //Q
    GranVolumeTable[6]=12.3 ; //E
    GranVolumeTable[7]=9.0  ; //G
    GranVolumeTable[8]=10.4 ; //H
    GranVolumeTable[9]=5.2  ; //I
    GranVolumeTable[10]=4.9 ; //L
    GranVolumeTable[11]=11.3; //K
    GranVolumeTable[12]=5.7 ; //M
    GranVolumeTable[13]=5.2 ; //F
    GranVolumeTable[14]=8.0 ; //P
    GranVolumeTable[15]=9.2 ; //S
    GranVolumeTable[16]=8.6 ; //T
    GranVolumeTable[17]=5.4 ; //W
    GranVolumeTable[18]=6.2 ; //Y
    GranVolumeTable[19]=5.9 ; //V
*/
}

MDOUBLE granthamChemicalDistances::getHughesHydrophobicityDistance(
		const int aa1,const int aa2) const {
	int v1=0;
	int v2=0;
	if ((aa1==0) || (aa1==4) || (aa1==13) || //acf
		(aa1==7) || (aa1==8) || (aa1==9) || //ghi
		(aa1==11) || (aa1==10) || (aa1==12) || //klm
		(aa1==16) || (aa1==19) || (aa1==17)
		 || (aa1==18))  //tvwy
		v1=1;
	if ((aa2==0) || (aa2==4) || (aa2==13) || //acf
		(aa2==7) || (aa2==8) || (aa2==9) || //ghi
		(aa2==11) || (aa2==10) || (aa2==12) || //klm
		(aa2==16) || (aa2==19) || (aa2==17)
		 || (aa2==18))  //tvwy
		v2=1;

	if (v1!=v2) return 1; 
	return 0;
}

MDOUBLE granthamChemicalDistances::getHughesPolarityDistance(
		const int aa1,const int aa2) const {
	int v1=0;
	int v2=0;
	if ((aa1==4) || (aa1==3) || (aa1==6) || //cde
		(aa1==8) || (aa1==11) || (aa1==2) || //hkn
		(aa1==5) || (aa1==1) || (aa1==15) || //qrs
		(aa1==16) || (aa1==17) || (aa1==18))  //tyw
		v1=1;
	if ((aa2==4) || (aa2==3) || (aa2==6) || //cde
		(aa2==8) || (aa2==11) || (aa2==2) || //hkn
		(aa2==5) || (aa2==1) || (aa2==15) || //qrs
		(aa2==16) || (aa2==17) || (aa2==18))  //tyw
		v2=1;

	if (v1!=v2) return 1; 
	return 0;
}
MDOUBLE granthamChemicalDistances::getHughesChargeDistance(
		const int aa1,const int aa2) const {
	int v1=0;
	int v2=0;
	if ((aa1==8) || (aa1==11) || (aa1==1)) v1=1;
	if ( (aa1==3) || (aa1==6)) v1=2;
	else v1=3;

	if ((aa2==8) || (aa2==11) || (aa2==1)) v2=1;
	if ( (aa2==3) || (aa2==6)) v2=2;
	else v2=3;

	if (v1!=v2) return 1; 
	return 0;
}



MDOUBLE granthamChemicalDistances::getGranthamDistance(const int aa1, const int aa2) const {
		if (aa1>aa2) return GranChemDist[aa2][aa1] ;
		else return GranChemDist[aa1][aa2];
}

MDOUBLE granthamChemicalDistances::getGranthamPolarityDistance(const int aa1,const int aa2) const{
	return fabs(GranPolarityTable[aa1]-GranPolarityTable[aa2]);
}

MDOUBLE granthamChemicalDistances::getGranthamPolarity(const int aa1) const{
	return GranPolarityTable[aa1];
}




