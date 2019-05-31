#include "PDE.h"
void SetTetraIntGaus4 ();
void SetTetraIntGaus15();
void SetTetraIntNode4 ();
void SetTetraIntNode10();

void SetPrismIntGaus6 ();
void SetPrismIntGaus18();
void SetPrismIntNode6 ();
void SetPrismIntNode18();

void SetHexahIntGaus8 ();
void SetHexahIntGaus27();
void SetHexahIntNode8 ();
void SetHexahIntNode20();
void SetHexahIntNode27();

void SetTrianIntGaus3();
void SetTrianIntGaus6();
void SetTrianIntNode3();
void SetTrianIntNode6();

void SetQuadrIntGaus4();
void SetQuadrIntGaus9();
void SetQuadrIntNode4();
void SetQuadrIntNode8();
void SetQuadrIntNode9();

void SetLineIntGaus2();
void SetLineIntGaus3();
void SetLineIntNode2();
void SetLineIntNode3();

void SetPointIntNode1();

void InitialGaus(int Type)
{
    switch(ElemType[Type])
    {
        case P1:
            RefcDim = 0;
            NodeNum = 1;
            SetPointIntNode1();
            break;
        case L2:
            RefcDim = 1;
            NodeNum = 2;
            switch(integtype)
            {
                case NodeInt:   SetLineIntNode2();    break;
                default:        SetLineIntGaus2();    break;
                case GausInt3:  SetLineIntGaus3();    break;
            }
            break;
        case L3:
            RefcDim = 1;
            NodeNum = 3;
            switch(integtype)
            {
                case NodeInt:   SetLineIntNode3();    break;
                default:        SetLineIntGaus3();    break;
            }
            break;
        case T3:
            RefcDim = 2;
            NodeNum = 3;
            switch(integtype)
            {
                case NodeInt:   SetTrianIntNode3();   break;
                default:        SetTrianIntGaus3();   break;
                case GausInt3:  SetTrianIntGaus6();   break;
            }
            break;
        case T6:
            RefcDim = 2;
            NodeNum = 6;
            switch(integtype)
            {
                case NodeInt:   SetTrianIntNode6();   break;
                default:        SetTrianIntGaus6();   break;
            }
            break;
        case Q4:
            RefcDim = 2;
            NodeNum = 4;
            switch(integtype)
            {
                case NodeInt:   SetQuadrIntNode4();   break;
                default:        SetQuadrIntGaus4();   break;
                case GausInt3:  SetQuadrIntGaus9();   break;
            }
            break;
        case Q8:
            RefcDim = 2;
            NodeNum = 8;
            switch(integtype)
            {
                case NodeInt:   SetQuadrIntNode8();   break;
                default:        SetQuadrIntGaus9();   break;
            }
            break;
        case Q9:
            RefcDim = 2;
            NodeNum = 9;
            switch(integtype)
            {
                case NodeInt:   SetQuadrIntNode9();   break;
                default:        SetQuadrIntGaus9();   break;
            }
            break;
        case W4:
            RefcDim = 3;
            NodeNum = 4;
            switch(integtype)
            {
                case NodeInt:   SetTetraIntNode4();   break;
                default:        SetTetraIntGaus4();   break;
                case GausInt3:  SetTetraIntGaus15();  break;
            }
            break;
        case W10:
            RefcDim = 3;
            NodeNum = 10;
            switch(integtype)
            {
                case NodeInt:   SetTetraIntNode10();  break;
                default:        SetTetraIntGaus15();  break;
            }
            break;
        case C8:
            RefcDim = 3;
            NodeNum = 8;
            switch(integtype)
            {
                case NodeInt:   SetHexahIntNode8();   break;
                default:        SetHexahIntGaus8();   break;
                case GausInt3:  SetHexahIntGaus27();  break;
            }
            break;
        case C20:
            RefcDim = 3;
            NodeNum = 20;
            switch(integtype)
            {
                case NodeInt:   SetHexahIntNode20();  break;
                default:        SetHexahIntGaus27();  break;
            }
            break;
        case C27:
            RefcDim = 3;
            NodeNum = 27;
            switch(integtype)
            {
                case NodeInt:   SetHexahIntNode27();  break;
                default:        SetHexahIntGaus27();  break;
            }
            break;
        case H6:
            RefcDim = 3;
            NodeNum = 6;
            switch(integtype)
            {
                case NodeInt:   SetPrismIntNode6();   break;
                default:        SetPrismIntGaus6();   break;
                case GausInt3:  SetPrismIntGaus18();  break;
            }
            break;
        case H15:
            break; // need write
            GausNum = 18;
            RefcDim = 3;
            NodeNum = 15;
            break;
        case H18:
            RefcDim = 3;
            NodeNum = 18;
            switch(integtype)
            {
                case NodeInt:   SetPrismIntNode18();  break;
                default:        SetPrismIntGaus18();  break;
            }
            break;
        case P5: ///
            GausNum = 1;
            RefcDim = 0;
            NodeNum = 1;
            break;
        case P13: ///
            GausNum = 1;
            RefcDim = 0;
            NodeNum = 1;
            break;
        case P14: ///
            GausNum = 1;
            RefcDim = 0;
            NodeNum = 1;
            break;
        default:
            break;
    }

}

void SetTetraIntGaus4()  //W
{
    GausNum = 4;
    
    GausCoor[(1-1)*GausNum + 1-1] = 0.5854101966;
    GausCoor[(2-1)*GausNum + 1-1] = 0.1381966011;
    GausCoor[(3-1)*GausNum + 1-1] = 0.1381966011;
    GausWeight[1-1] = 1.0/24;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.1381966011;
    GausCoor[(2-1)*GausNum + 2-1] = 0.5854101966;
    GausCoor[(3-1)*GausNum + 2-1] = 0.1381966011;
    GausWeight[2-1] = 1.0/24;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.1381966011;
    GausCoor[(2-1)*GausNum + 3-1] = 0.1381966011;
    GausCoor[(3-1)*GausNum + 3-1] = 0.5854101966;
    GausWeight[3-1] = 1.0/24;
    
    GausCoor[(1-1)*GausNum + 4-1] = 0.1381966011;
    GausCoor[(2-1)*GausNum + 4-1] = 0.1381966011;
    GausCoor[(3-1)*GausNum + 4-1] = 0.1381966011;
    GausWeight[4-1] = 1.0/24;
    
    return;
}
void SetTetraIntGaus15() //W
{
    GausNum = 15;
    
    GausCoor[(1-1)*GausNum +  1-1] = 0.2500000000;
    GausCoor[(2-1)*GausNum +  1-1] = 0.2500000000;
    GausCoor[(3-1)*GausNum +  1-1] = 0.2500000000;
    GausWeight[1-1] = 0.0197530864;
    
    GausCoor[(1-1)*GausNum +  2-1] = 0.0919710780;
    GausCoor[(2-1)*GausNum +  2-1] = 0.0919710780;
    GausCoor[(3-1)*GausNum +  2-1] = 0.0919710780;
    GausWeight[2-1] = 0.0119895139;
    
    GausCoor[(1-1)*GausNum +  3-1] = 0.7240867658;
    GausCoor[(2-1)*GausNum +  3-1] = 0.0919710780;
    GausCoor[(3-1)*GausNum +  3-1] = 0.0919710780;
    GausWeight[3-1] = 0.0119895139;
    
    GausCoor[(1-1)*GausNum +  4-1] = 0.0919710780;
    GausCoor[(2-1)*GausNum +  4-1] = 0.7240867658;
    GausCoor[(3-1)*GausNum +  4-1] = 0.0919710780;
    GausWeight[4-1] = 0.0119895139;
    
    GausCoor[(1-1)*GausNum +  5-1] = 0.0919710780;
    GausCoor[(2-1)*GausNum +  5-1] = 0.0919710780;
    GausCoor[(3-1)*GausNum +  5-1] = 0.7240867658;
    GausWeight[5-1] = 0.0119895139;
    
    GausCoor[(1-1)*GausNum +  6-1] = 0.3197936278;
    GausCoor[(2-1)*GausNum +  6-1] = 0.3197936278;
    GausCoor[(3-1)*GausNum +  6-1] = 0.3197936278;
    GausWeight[6-1] = 0.0115113678;
    
    GausCoor[(1-1)*GausNum +  7-1] = 0.0406191165;
    GausCoor[(2-1)*GausNum +  7-1] = 0.3197936278;
    GausCoor[(3-1)*GausNum +  7-1] = 0.3197936278;
    GausWeight[7-1] = 0.0115113678;
    
    GausCoor[(1-1)*GausNum +  8-1] = 0.3197936278;
    GausCoor[(2-1)*GausNum +  8-1] = 0.0406191165;
    GausCoor[(3-1)*GausNum +  8-1] = 0.3197936278;
    GausWeight[8-1] = 0.0115113678;
    
    GausCoor[(1-1)*GausNum +  9-1] = 0.3197936278;
    GausCoor[(2-1)*GausNum +  9-1] = 0.3197936278;
    GausCoor[(3-1)*GausNum +  9-1] = 0.0406191165;
    GausWeight[9-1] = 0.0115113678;
    
    GausCoor[(1-1)*GausNum + 10-1] = 0.4436491673;
    GausCoor[(2-1)*GausNum + 10-1] = 0.0563508326;
    GausCoor[(3-1)*GausNum + 10-1] = 0.0563508326;
    GausWeight[10-1] = 0.0088183421;
    
    GausCoor[(1-1)*GausNum + 11-1] = 0.0563508326;
    GausCoor[(2-1)*GausNum + 11-1] = 0.4436491673;
    GausCoor[(3-1)*GausNum + 11-1] = 0.0563508326;
    GausWeight[11-1] = 0.0088183421;
    
    GausCoor[(1-1)*GausNum + 12-1] = 0.0563508326;
    GausCoor[(2-1)*GausNum + 12-1] = 0.0563508326;
    GausCoor[(3-1)*GausNum + 12-1] = 0.4436491673;
    GausWeight[12-1] = 0.0088183421;
    
    GausCoor[(1-1)*GausNum + 13-1] = 0.4436491673;
    GausCoor[(2-1)*GausNum + 13-1] = 0.0563508326;
    GausCoor[(3-1)*GausNum + 13-1] = 0.4436491673;
    GausWeight[13-1] = 0.0088183421;
    
    GausCoor[(1-1)*GausNum + 14-1] = 0.4436491673;
    GausCoor[(2-1)*GausNum + 14-1] = 0.4436491673;
    GausCoor[(3-1)*GausNum + 14-1] = 0.0563508326;
    GausWeight[14-1] = 0.0088183421;
    
    GausCoor[(1-1)*GausNum + 15-1] = 0.0563508326;
    GausCoor[(2-1)*GausNum + 15-1] = 0.4436491673;
    GausCoor[(3-1)*GausNum + 15-1] = 0.4436491673;
    GausWeight[15-1] = 0.0088183421;
    
    return;
}
void SetTetraIntNode4()  //W
{
    GausNum = 4;
    
    GausCoor[(1-1)*GausNum + 1-1] = 1.0;
    GausCoor[(2-1)*GausNum + 1-1] = 0.0;
    GausCoor[(3-1)*GausNum + 1-1] = 0.0;
    GausWeight[1-1] = 0.041666667;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.0;
    GausCoor[(2-1)*GausNum + 2-1] = 1.0;
    GausCoor[(3-1)*GausNum + 2-1] = 0.0;
    GausWeight[2-1] = 0.041666667;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.0;
    GausCoor[(2-1)*GausNum + 3-1] = 0.0;
    GausCoor[(3-1)*GausNum + 3-1] = 0.0;
    GausWeight[3-1] = 0.041666667;
    
    GausCoor[(1-1)*GausNum + 4-1] = 0.0;
    GausCoor[(2-1)*GausNum + 4-1] = 0.0;
    GausCoor[(3-1)*GausNum + 4-1] = 1.0;
    GausWeight[4-1] = 0.041666667;
    
    return;
}
void SetTetraIntNode10() //W
{
    GausNum = 10;
    
    GausCoor[(1-1)*GausNum +  1-1] = 1.0;
    GausCoor[(2-1)*GausNum +  1-1] = 0.0;
    GausCoor[(3-1)*GausNum +  1-1] = 0.0;
    GausWeight[1-1] = 1.0/96;
    
    GausCoor[(1-1)*GausNum +  2-1] = 0.5;
    GausCoor[(2-1)*GausNum +  2-1] = 0.5;
    GausCoor[(3-1)*GausNum +  2-1] = 0.0;
    GausWeight[2-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum +  3-1] = 0.0;
    GausCoor[(2-1)*GausNum +  3-1] = 1.0;
    GausCoor[(3-1)*GausNum +  3-1] = 0.0;
    GausWeight[3-1] = 1.0/96;
    
    GausCoor[(1-1)*GausNum +  4-1] = 0.0;
    GausCoor[(2-1)*GausNum +  4-1] = 0.5;
    GausCoor[(3-1)*GausNum +  4-1] = 0.5;
    GausWeight[4-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum +  5-1] = 0.0;
    GausCoor[(2-1)*GausNum +  5-1] = 0.0;
    GausCoor[(3-1)*GausNum +  5-1] = 1.0;
    GausWeight[5-1] = 1.0/96;
    
    GausCoor[(1-1)*GausNum +  6-1] = 0.5;
    GausCoor[(2-1)*GausNum +  6-1] = 0.0;
    GausCoor[(3-1)*GausNum +  6-1] = 0.5;
    GausWeight[6-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum +  7-1] = 0.5;
    GausCoor[(2-1)*GausNum +  7-1] = 0.0;
    GausCoor[(3-1)*GausNum +  7-1] = 0.0;
    GausWeight[7-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum +  8-1] = 0.0;
    GausCoor[(2-1)*GausNum +  8-1] = 0.5;
    GausCoor[(3-1)*GausNum +  8-1] = 0.0;
    GausWeight[8-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum +  9-1] = 0.0;
    GausCoor[(2-1)*GausNum +  9-1] = 0.0;
    GausCoor[(3-1)*GausNum +  9-1] = 0.5;
    GausWeight[9-1] = 1.0/48;
    
    GausCoor[(1-1)*GausNum + 10-1] = 0.0;
    GausCoor[(2-1)*GausNum + 10-1] = 0.0;
    GausCoor[(3-1)*GausNum + 10-1] = 0.0;
    GausWeight[10-1] = 1.0/96;
    
    return;
}
void SetPrismIntGaus6()  //H
{
    GausNum = 6;
    
    GausCoor[(1-1)*GausNum + 1-1] = 2.0/3.0;
    GausCoor[(2-1)*GausNum + 1-1] = 1.0/6.0;
    GausCoor[(3-1)*GausNum + 1-1] = 0.5773502691;
    GausWeight[1-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 2-1] = 2.0/3.0;
    GausCoor[(2-1)*GausNum + 2-1] = 1.0/6.0;
    GausCoor[(3-1)*GausNum + 2-1] = -0.577350269;
    GausWeight[2-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 3-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 3-1] = 1.0/6.0;
    GausCoor[(3-1)*GausNum + 3-1] = 0.5773502691;
    GausWeight[3-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 4-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 4-1] = 1.0/6.0;
    GausCoor[(3-1)*GausNum + 4-1] = -0.577350269;
    GausWeight[4-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 5-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 5-1] = 2.0/3.0;
    GausCoor[(3-1)*GausNum + 5-1] = 0.5773502691;
    GausWeight[5-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 6-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 6-1] = 2.0/3.0;
    GausCoor[(3-1)*GausNum + 6-1] = -0.577350269;
    GausWeight[6-1] = 1.0/6.0;
    
    return;
}
void SetPrismIntGaus18() //H
{
    GausNum = 18;
    
    GausCoor[(1-1)*GausNum + 1-1] = 0.8168475729;
    GausCoor[(2-1)*GausNum + 1-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 1-1] = 0.7745966692;
    GausWeight[1-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.8168475729;
    GausCoor[(2-1)*GausNum + 2-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 2-1] = 0.0000000000;
    GausWeight[2-1] = 0.8888888888*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.8168475729;
    GausCoor[(2-1)*GausNum + 3-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 3-1] = -0.774596669;
    GausWeight[3-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 4-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 4-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 4-1] = 0.7745966692;
    GausWeight[4-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 5-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 5-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 5-1] = 0.0000000000;
    GausWeight[5-1] = 0.8888888888*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 6-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 6-1] = 0.0915762135;
    GausCoor[(3-1)*GausNum + 6-1] = -0.774596669;
    GausWeight[6-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 7-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 7-1] = 0.8168475729;
    GausCoor[(3-1)*GausNum + 7-1] = 0.7745966692;
    GausWeight[7-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 8-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 8-1] = 0.8168475729;
    GausCoor[(3-1)*GausNum + 8-1] = 0.0000000000;
    GausWeight[8-1] = 0.8888888888*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 9-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 9-1] = 0.8168475729;
    GausCoor[(3-1)*GausNum + 9-1] = -0.774596669;
    GausWeight[9-1] = 0.5555555555*0.0549758718;
    
    GausCoor[(1-1)*GausNum + 10-1] = 0.1081030181;
    GausCoor[(2-1)*GausNum + 10-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 10-1] = 0.7745966692;
    GausWeight[10-1] = 0.5555555555*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 11-1] = 0.1081030181;
    GausCoor[(2-1)*GausNum + 11-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 11-1] = 0.0000000000;
    GausWeight[11-1] = 0.8888888888*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 12-1] = 0.1081030181;
    GausCoor[(2-1)*GausNum + 12-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 12-1] = -0.774596669;
    GausWeight[12-1] = 0.5555555555*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 13-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 13-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 13-1] = 0.7745966692;
    GausWeight[13-1] = 0.5555555555*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 14-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 14-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 14-1] = 0.0000000000;
    GausWeight[14-1] = 0.8888888888*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 15-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 15-1] = 0.4459484909;
    GausCoor[(3-1)*GausNum + 15-1] = -0.774596669;
    GausWeight[15-1] = 0.5555555555*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 16-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 16-1] = 0.1081030181;
    GausCoor[(3-1)*GausNum + 16-1] = 0.7745966692;
    GausWeight[16-1] = 0.5555555555*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 17-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 17-1] = 0.1081030181;
    GausCoor[(3-1)*GausNum + 17-1] = 0.0000000000;
    GausWeight[17-1] = 0.8888888888*0.1116907948;
    
    GausCoor[(1-1)*GausNum + 18-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 18-1] = 0.1081030181;
    GausCoor[(3-1)*GausNum + 18-1] = -0.774596669;
    GausWeight[18-1] = 0.5555555555*0.1116907948;
    
    return;
}
void SetPrismIntNode6()  //H
{
    GausNum = 6;
    
    GausCoor[(1-1)*GausNum + 1-1] =  1.0;
    GausCoor[(2-1)*GausNum + 1-1] =  0.0;
    GausCoor[(3-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 1.666666e-01;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.0;
    GausCoor[(2-1)*GausNum + 2-1] =  1.0;
    GausCoor[(3-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 1.666666e-01;
    
    GausCoor[(1-1)*GausNum + 3-1] =  0.0;
    GausCoor[(2-1)*GausNum + 3-1] =  0.0;
    GausCoor[(3-1)*GausNum + 3-1] = -1.0;
    GausWeight[3-1] = 1.666666e-01;
    
    GausCoor[(1-1)*GausNum + 4-1] =  1.0;
    GausCoor[(2-1)*GausNum + 4-1] =  0.0;
    GausCoor[(3-1)*GausNum + 4-1] =  1.0;
    GausWeight[4-1] = 1.666666e-01;
    
    GausCoor[(1-1)*GausNum + 5-1] =  0.0;
    GausCoor[(2-1)*GausNum + 5-1] =  1.0;
    GausCoor[(3-1)*GausNum + 5-1] =  1.0;
    GausWeight[5-1] = 1.666666e-01;
    
    GausCoor[(1-1)*GausNum + 6-1] =  0.0;
    GausCoor[(2-1)*GausNum + 6-1] =  0.0;
    GausCoor[(3-1)*GausNum + 6-1] =  1.0;
    GausWeight[6-1] = 1.666666e-01;
    
    return;
}
void SetPrismIntNode18() //H
{
    GausNum = 18;
    
    GausCoor[(1-1)*GausNum + 1-1] =  1.0;
    GausCoor[(2-1)*GausNum + 1-1] =  0.0;
    GausCoor[(3-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.5;
    GausCoor[(2-1)*GausNum + 2-1] =  0.5;
    GausCoor[(3-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 6.250000e-02;
    
    GausCoor[(1-1)*GausNum + 3-1] =  0.0;
    GausCoor[(2-1)*GausNum + 3-1] =  1.0;
    GausCoor[(3-1)*GausNum + 3-1] = -1.0;
    GausWeight[3-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 4-1] =  0.0;
    GausCoor[(2-1)*GausNum + 4-1] =  0.5;
    GausCoor[(3-1)*GausNum + 4-1] = -1.0;
    GausWeight[4-1] = 6.250000e-02;
    
    GausCoor[(1-1)*GausNum + 5-1] =  0.0;
    GausCoor[(2-1)*GausNum + 5-1] =  0.0;
    GausCoor[(3-1)*GausNum + 5-1] = -1.0;
    GausWeight[5-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 6-1] =  0.5;
    GausCoor[(2-1)*GausNum + 6-1] =  0.0;
    GausCoor[(3-1)*GausNum + 6-1] = -1.0;
    GausWeight[6-1] = 6.250000e-02;
    
    GausCoor[(1-1)*GausNum + 7-1] =  1.0;
    GausCoor[(2-1)*GausNum + 7-1] =  0.0;
    GausCoor[(3-1)*GausNum + 7-1] =  0.0;
    GausWeight[7-1] = 4.166667e-02;
    
    GausCoor[(1-1)*GausNum + 8-1] =  0.5;
    GausCoor[(2-1)*GausNum + 8-1] =  0.5;
    GausCoor[(3-1)*GausNum + 8-1] =  0.0;
    GausWeight[8-1] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 9-1] =  0.0;
    GausCoor[(2-1)*GausNum + 9-1] =  1.0;
    GausCoor[(3-1)*GausNum + 9-1] =  0.0;
    GausWeight[9-1] = 4.166667e-02;
    
    GausCoor[(1-1)*GausNum + 10-1] = 0.0;
    GausCoor[(2-1)*GausNum + 10-1] = 0.5;
    GausCoor[(3-1)*GausNum + 10-1] = 0.0;
    GausWeight[10-1] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 11-1] = 0.0;
    GausCoor[(2-1)*GausNum + 11-1] = 0.0;
    GausCoor[(3-1)*GausNum + 11-1] = 0.0;
    GausWeight[11-1] = 4.166667e-02;
    
    GausCoor[(1-1)*GausNum + 12-1] = 0.5;
    GausCoor[(2-1)*GausNum + 12-1] = 0.0;
    GausCoor[(3-1)*GausNum + 12-1] = 0.0;
    GausWeight[12-1] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 13-1] = 1.0;
    GausCoor[(2-1)*GausNum + 13-1] = 0.0;
    GausCoor[(3-1)*GausNum + 13-1] = 1.0;
    GausWeight[13-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 14-1] = 0.5;
    GausCoor[(2-1)*GausNum + 14-1] = 0.5;
    GausCoor[(3-1)*GausNum + 14-1] = 1.0;
    GausWeight[14-1] = 6.250000e-02;
    
    GausCoor[(1-1)*GausNum + 15-1] = 0.0;
    GausCoor[(2-1)*GausNum + 15-1] = 1.0;
    GausCoor[(3-1)*GausNum + 15-1] = 1.0;
    GausWeight[15-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 16-1] = 0.0;
    GausCoor[(2-1)*GausNum + 16-1] = 0.5;
    GausCoor[(3-1)*GausNum + 16-1] = 1.0;
    GausWeight[16-1] = 6.250000e-02;
    
    GausCoor[(1-1)*GausNum + 17-1] = 0.0;
    GausCoor[(2-1)*GausNum + 17-1] = 0.0;
    GausCoor[(3-1)*GausNum + 17-1] = 1.0;
    GausWeight[17-1] = 2.083334e-02;
    
    GausCoor[(1-1)*GausNum + 18-1] = 0.5;
    GausCoor[(2-1)*GausNum + 18-1] = 0.0;
    GausCoor[(3-1)*GausNum + 18-1] = 1.0;
    GausWeight[18-1] = 6.250000e-02;
    
    return;
}
void SetHexahIntGaus8()  //C
{
    GausNum = 8;
    
    GausCoor[(1-1)*GausNum + 1-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 1-1] =  5.773502692e-001;
    GausCoor[(3-1)*GausNum + 1-1] =  5.773502692e-001;
    GausWeight[1-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 2-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 2-1] =  5.773502692e-001;
    GausCoor[(3-1)*GausNum + 2-1] = -5.773502692e-001;
    GausWeight[2-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 3-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 3-1] = -5.773502692e-001;
    GausCoor[(3-1)*GausNum + 3-1] =  5.773502692e-001;
    GausWeight[3-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 4-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 4-1] = -5.773502692e-001;
    GausCoor[(3-1)*GausNum + 4-1] = -5.773502692e-001;
    GausWeight[4-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 5-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 5-1] =  5.773502692e-001;
    GausCoor[(3-1)*GausNum + 5-1] =  5.773502692e-001;
    GausWeight[5-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 6-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 6-1] =  5.773502692e-001;
    GausCoor[(3-1)*GausNum + 6-1] = -5.773502692e-001;
    GausWeight[6-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 7-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 7-1] = -5.773502692e-001;
    GausCoor[(3-1)*GausNum + 7-1] =  5.773502692e-001;
    GausWeight[7-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 8-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 8-1] = -5.773502692e-001;
    GausCoor[(3-1)*GausNum + 8-1] = -5.773502692e-001;
    GausWeight[8-1] = 1.0;
    
    return;
}
void SetHexahIntGaus27() //C
{
    GausNum = 27;
    
    GausCoor[(1-1)*GausNum +  1-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  1-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum +  1-1] =  7.745966692e-001;
    GausWeight[1-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum +  2-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  2-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum +  2-1] =  0.000000000e+000;
    GausWeight[2-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum +  3-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  3-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum +  3-1] = -7.745966692e-001;
    GausWeight[3-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum +  4-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  4-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum +  4-1] =  7.745966692e-001;
    GausWeight[4-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum +  5-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  5-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum +  5-1] =  0.000000000e+000;
    GausWeight[5-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum +  6-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  6-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum +  6-1] = -7.745966692e-001;
    GausWeight[6-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum +  7-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  7-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum +  7-1] =  7.745966692e-001;
    GausWeight[7-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum +  8-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  8-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum +  8-1] =  0.000000000e+000;
    GausWeight[8-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum +  9-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum +  9-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum +  9-1] = -7.745966692e-001;
    GausWeight[9-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum + 10-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 10-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 10-1] =  7.745966692e-001;
    GausWeight[10-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 11-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 11-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 11-1] =  0.000000000e+000;
    GausWeight[11-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum + 12-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 12-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 12-1] = -7.745966692e-001;
    GausWeight[12-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 13-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 13-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 13-1] =  7.745966692e-001;
    GausWeight[13-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum + 14-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 14-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 14-1] =  0.000000000e+000;
    GausWeight[14-1] = 7.023319616e-001;
    
    GausCoor[(1-1)*GausNum + 15-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 15-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 15-1] = -7.745966692e-001;
    GausWeight[15-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum + 16-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 16-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 16-1] =  7.745966692e-001;
    GausWeight[16-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 17-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 17-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 17-1] =  0.000000000e+000;
    GausWeight[17-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum + 18-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 18-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 18-1] = -7.745966692e-001;
    GausWeight[18-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 19-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 19-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 19-1] =  7.745966692e-001;
    GausWeight[19-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum + 20-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 20-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 20-1] =  0.000000000e+000;
    GausWeight[20-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 21-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 21-1] =  7.745966692e-001;
    GausCoor[(3-1)*GausNum + 21-1] = -7.745966692e-001;
    GausWeight[21-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum + 22-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 22-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 22-1] =  7.745966692e-001;
    GausWeight[22-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 23-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 23-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 23-1] =  0.000000000e+000;
    GausWeight[23-1] = 4.389574760e-001;
    
    GausCoor[(1-1)*GausNum + 24-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 24-1] =  0.000000000e+000;
    GausCoor[(3-1)*GausNum + 24-1] = -7.745966692e-001;
    GausWeight[24-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 25-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 25-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 25-1] =  7.745966692e-001;
    GausWeight[25-1] = 1.714677641e-001;
    
    GausCoor[(1-1)*GausNum + 26-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 26-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 26-1] =  0.000000000e+000;
    GausWeight[26-1] = 2.743484225e-001;
    
    GausCoor[(1-1)*GausNum + 27-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 27-1] = -7.745966692e-001;
    GausCoor[(3-1)*GausNum + 27-1] = -7.745966692e-001;
    GausWeight[27-1] = 1.714677641e-001;
    
    return;
}
void SetHexahIntNode8()  //C
{
    GausNum = 8;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausCoor[(2-1)*GausNum + 1-1] = -1.0;
    GausCoor[(3-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 2-1] =  1.0;
    GausCoor[(2-1)*GausNum + 2-1] = -1.0;
    GausCoor[(3-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 3-1] =  1.0;
    GausCoor[(2-1)*GausNum + 3-1] =  1.0;
    GausCoor[(3-1)*GausNum + 3-1] = -1.0;
    GausWeight[3-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 4-1] = -1.0;
    GausCoor[(2-1)*GausNum + 4-1] =  1.0;
    GausCoor[(3-1)*GausNum + 4-1] = -1.0;
    GausWeight[4-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 5-1] = -1.0;
    GausCoor[(2-1)*GausNum + 5-1] = -1.0;
    GausCoor[(3-1)*GausNum + 5-1] =  1.0;
    GausWeight[5-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 6-1] =  1.0;
    GausCoor[(2-1)*GausNum + 6-1] = -1.0;
    GausCoor[(3-1)*GausNum + 6-1] =  1.0;
    GausWeight[6-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 7-1] =  1.0;
    GausCoor[(2-1)*GausNum + 7-1] =  1.0;
    GausCoor[(3-1)*GausNum + 7-1] =  1.0;
    GausWeight[7-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 8-1] = -1.0;
    GausCoor[(2-1)*GausNum + 8-1] =  1.0;
    GausCoor[(3-1)*GausNum + 8-1] =  1.0;
    GausWeight[8-1] = 1.0;
    
    return;
}
void SetHexahIntNode20() //C
{
    GausNum = 20;
    
    GausCoor[(1-1)*GausNum +  1-1] = -1.0;
    GausCoor[(2-1)*GausNum +  1-1] = -1.0;
    GausCoor[(3-1)*GausNum +  1-1] = -1.0;
    GausWeight[1] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  2-1] =  0.0;
    GausCoor[(2-1)*GausNum +  2-1] = -1.0;
    GausCoor[(3-1)*GausNum +  2-1] = -1.0;
    GausWeight[2] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum +  3-1] =  1.0;
    GausCoor[(2-1)*GausNum +  3-1] = -1.0;
    GausCoor[(3-1)*GausNum +  3-1] = -1.0;
    GausWeight[3] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  4-1] = -1.0;
    GausCoor[(2-1)*GausNum +  4-1] =  0.0;
    GausCoor[(3-1)*GausNum +  4-1] = -1.0;
    GausWeight[4] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum +  5-1] =  1.0;
    GausCoor[(2-1)*GausNum +  5-1] =  0.0;
    GausCoor[(3-1)*GausNum +  5-1] = -1.0;
    GausWeight[5] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum +  6-1] = -1.0;
    GausCoor[(2-1)*GausNum +  6-1] =  1.0;
    GausCoor[(3-1)*GausNum +  6-1] = -1.0;
    GausWeight[6] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  7-1] =  0.0;
    GausCoor[(2-1)*GausNum +  7-1] =  1.0;
    GausCoor[(3-1)*GausNum +  7-1] = -1.0;
    GausWeight[7] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum +  8-1] =  1.0;
    GausCoor[(2-1)*GausNum +  8-1] =  1.0;
    GausCoor[(3-1)*GausNum +  8-1] = -1.0;
    GausWeight[8] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  9-1] = -1.0;
    GausCoor[(2-1)*GausNum +  9-1] = -1.0;
    GausCoor[(3-1)*GausNum +  9-1] =  0.0;
    GausWeight[9] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 10-1] =  1.0;
    GausCoor[(2-1)*GausNum + 10-1] = -1.0;
    GausCoor[(3-1)*GausNum + 10-1] =  0.0;
    GausWeight[10] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 11-1] = -1.0;
    GausCoor[(2-1)*GausNum + 11-1] =  1.0;
    GausCoor[(3-1)*GausNum + 11-1] =  0.0;
    GausWeight[11] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 12-1] =  1.0;
    GausCoor[(2-1)*GausNum + 12-1] =  1.0;
    GausCoor[(3-1)*GausNum + 12-1] =  0.0;
    GausWeight[12] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 13-1] = -1.0;
    GausCoor[(2-1)*GausNum + 13-1] = -1.0;
    GausCoor[(3-1)*GausNum + 13-1] =  1.0;
    GausWeight[13] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 14-1] =  0.0;
    GausCoor[(2-1)*GausNum + 14-1] = -1.0;
    GausCoor[(3-1)*GausNum + 14-1] =  1.0;
    GausWeight[14] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 15-1] =  1.0;
    GausCoor[(2-1)*GausNum + 15-1] = -1.0;
    GausCoor[(3-1)*GausNum + 15-1] =  1.0;
    GausWeight[15] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 16-1] = -1.0;
    GausCoor[(2-1)*GausNum + 16-1] =  0.0;
    GausCoor[(3-1)*GausNum + 16-1] =  1.0;
    GausWeight[16] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 17-1] =  1.0;
    GausCoor[(2-1)*GausNum + 17-1] =  0.0;
    GausCoor[(3-1)*GausNum + 17-1] =  1.0;
    GausWeight[17] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 18-1] = -1.0;
    GausCoor[(2-1)*GausNum + 18-1] =  1.0;
    GausCoor[(3-1)*GausNum + 18-1] =  1.0;
    GausWeight[18] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 19-1] =  0.0;
    GausCoor[(2-1)*GausNum + 19-1] =  1.0;
    GausCoor[(3-1)*GausNum + 19-1] =  1.0;
    GausWeight[19] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 20-1] =  1.0;
    GausCoor[(2-1)*GausNum + 20-1] =  1.0;
    GausCoor[(3-1)*GausNum + 20-1] =  1.0;
    GausWeight[20] = 2.500000e-01;
    
    return;
}
void SetHexahIntNode27() //C
{
    GausNum = 27;
    
    GausCoor[(1-1)*GausNum +  1-1] = -1.0;
    GausCoor[(2-1)*GausNum +  1-1] = -1.0;
    GausCoor[(3-1)*GausNum +  1-1] = -1.0;
    GausWeight[1] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum +  2-1] =  0.0;
    GausCoor[(2-1)*GausNum +  2-1] = -1.0;
    GausCoor[(3-1)*GausNum +  2-1] = -1.0;
    GausWeight[2] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  3-1] =  1.0;
    GausCoor[(2-1)*GausNum +  3-1] = -1.0;
    GausCoor[(3-1)*GausNum +  3-1] = -1.0;
    GausWeight[3] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum +  4-1] = -1.0;
    GausCoor[(2-1)*GausNum +  4-1] =  0.0;
    GausCoor[(3-1)*GausNum +  4-1] = -1.0;
    GausWeight[4] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  5-1] =  0.0;
    GausCoor[(2-1)*GausNum +  5-1] =  0.0;
    GausCoor[(3-1)*GausNum +  5-1] = -1.0;
    GausWeight[5] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum +  6-1] =  1.0;
    GausCoor[(2-1)*GausNum +  6-1] =  0.0;
    GausCoor[(3-1)*GausNum +  6-1] = -1.0;
    GausWeight[6] = 2.500000e-01; 
    
    GausCoor[(1-1)*GausNum +  7-1] = -1.0;
    GausCoor[(2-1)*GausNum +  7-1] =  1.0;
    GausCoor[(3-1)*GausNum +  7-1] = -1.0;
    GausWeight[7] = 1.250000e-01;
    GausCoor[(1-1)*GausNum +  8-1] =  0.0;
    GausCoor[(2-1)*GausNum +  8-1] =  1.0;
    GausCoor[(3-1)*GausNum +  8-1] = -1.0;
    GausWeight[8] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum +  9-1] =  1.0;
    GausCoor[(2-1)*GausNum +  9-1] =  1.0;
    GausCoor[(3-1)*GausNum +  9-1] = -1.0;
    GausWeight[9] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 10-1] = -1.0;
    GausCoor[(2-1)*GausNum + 10-1] = -1.0;
    GausCoor[(3-1)*GausNum + 10-1] =  0.0;
    GausWeight[10] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 11-1] =  0.0;
    GausCoor[(2-1)*GausNum + 11-1] = -1.0;
    GausCoor[(3-1)*GausNum + 11-1] =  0.0;
    GausWeight[11] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 12-1] =  1.0;
    GausCoor[(2-1)*GausNum + 12-1] = -1.0;
    GausCoor[(3-1)*GausNum + 12-1] =  0.0;
    GausWeight[12] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 13-1] = -1.0;
    GausCoor[(2-1)*GausNum + 13-1] =  0.0;
    GausCoor[(3-1)*GausNum + 13-1] =  0.0;
    GausWeight[13] = 5.000000e-01;
    GausCoor[(1-1)*GausNum + 14-1] =  0.0;
    GausCoor[(2-1)*GausNum + 14-1] =  0.0;
    GausCoor[(3-1)*GausNum + 14-1] =  0.0;
    GausWeight[14] = 1.0;
    
    GausCoor[(1-1)*GausNum + 15-1] =  1.0;
    GausCoor[(2-1)*GausNum + 15-1] =  0.0;
    GausCoor[(3-1)*GausNum + 15-1] =  0.0;
    GausWeight[15] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 16-1] = -1.0;
    GausCoor[(2-1)*GausNum + 16-1] =  1.0;
    GausCoor[(3-1)*GausNum + 16-1] =  0.0;
    GausWeight[16] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 17-1] =  0.0;
    GausCoor[(2-1)*GausNum + 17-1] =  1.0;
    GausCoor[(3-1)*GausNum + 17-1] =  0.0;
    GausWeight[17] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 18-1] =  1.0;
    GausCoor[(2-1)*GausNum + 18-1] =  1.0;
    GausCoor[(3-1)*GausNum + 18-1] =  0.0;
    GausWeight[18] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 19-1] = -1.0;
    GausCoor[(2-1)*GausNum + 19-1] = -1.0;
    GausCoor[(3-1)*GausNum + 19-1] =  1.0;
    GausWeight[19] = 1.250000e-01;
    GausCoor[(1-1)*GausNum + 20-1] =  0.0;
    GausCoor[(2-1)*GausNum + 20-1] = -1.0;
    GausCoor[(3-1)*GausNum + 20-1] =  1.0;
    GausWeight[20] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 21-1] =  1.0;
    GausCoor[(2-1)*GausNum + 21-1] = -1.0;
    GausCoor[(3-1)*GausNum + 21-1] =  1.0;
    GausWeight[21] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 22-1] = -1.0;
    GausCoor[(2-1)*GausNum + 22-1] =  0.0;
    GausCoor[(3-1)*GausNum + 22-1] =  1.0;
    GausWeight[22] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 23-1] =  0.0;
    GausCoor[(2-1)*GausNum + 23-1] =  0.0;
    GausCoor[(3-1)*GausNum + 23-1] =  1.0;
    GausWeight[23] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 24-1] =  1.0;
    GausCoor[(2-1)*GausNum + 24-1] =  0.0;
    GausCoor[(3-1)*GausNum + 24-1] =  1.0;
    GausWeight[24] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 25-1] = -1.0;
    GausCoor[(2-1)*GausNum + 25-1] =  1.0;
    GausCoor[(3-1)*GausNum + 25-1] =  1.0;
    GausWeight[25] = 1.250000e-01;
    
    GausCoor[(1-1)*GausNum + 26-1] =  0.0;
    GausCoor[(2-1)*GausNum + 26-1] =  1.0;
    GausCoor[(3-1)*GausNum + 26-1] =  1.0;
    GausWeight[26] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 27-1] =  1.0;
    GausCoor[(2-1)*GausNum + 27-1] =  1.0;
    GausCoor[(3-1)*GausNum + 27-1] =  1.0;
    GausWeight[27] = 1.250000e-01;
    
    return;
}
void SetTrianIntGaus3()  //T
{
    GausNum = 3;
    
    GausCoor[(1-1)*GausNum + 1-1] = 2.0/3.0;
    GausCoor[(2-1)*GausNum + 1-1] = 1.0/6.0;
    GausWeight[1-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 2-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 2-1] = 1.0/6.0;
    GausWeight[2-1] = 1.0/6.0;
    
    GausCoor[(1-1)*GausNum + 3-1] = 1.0/6.0;
    GausCoor[(2-1)*GausNum + 3-1] = 2.0/3.0;
    GausWeight[3-1] = 1.0/6.0;
    
    return;
}
void SetTrianIntGaus6()  //T
{
    GausNum = 6;
    
    GausCoor[(1-1)*GausNum + 1-1] = 0.8168475729;
    GausCoor[(2-1)*GausNum + 1-1] = 0.0915762135;
    GausWeight[1-1] = 0.0549758718;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 2-1] = 0.0915762135;
    GausWeight[2-1] = 0.0549758718;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.0915762135;
    GausCoor[(2-1)*GausNum + 3-1] = 0.8168475729;
    GausWeight[3-1] = 0.0549758718;
    
    GausCoor[(1-1)*GausNum + 4-1] = 0.1081030181;
    GausCoor[(2-1)*GausNum + 4-1] = 0.4459484909;
    GausWeight[4-1] = 0.1116907948;
    
    GausCoor[(1-1)*GausNum + 5-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 5-1] = 0.4459484909;
    GausWeight[5-1] = 0.1116907948;
    
    GausCoor[(1-1)*GausNum + 6-1] = 0.4459484909;
    GausCoor[(2-1)*GausNum + 6-1] = 0.1081030181;
    GausWeight[6-1] = 0.1116907948;
    
    return;
}
void SetTrianIntNode3()  //T
{
    GausNum = 3;
    
    GausCoor[(1-1)*GausNum + 1-1] = 1.0;
    GausCoor[(2-1)*GausNum + 1-1] = 0.0;
    GausWeight[1-1] = 0.1666666;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.0;
    GausCoor[(2-1)*GausNum + 2-1] = 1.0;
    GausWeight[2-1] = 0.1666666;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.0;
    GausCoor[(2-1)*GausNum + 3-1] = 0.0;
    GausWeight[3-1] = 0.1666666;
    
    return;
}
void SetTrianIntNode6()  //T
{
    GausNum = 6;
    
    GausCoor[(1-1)*GausNum + 1-1] = 1.0;
    GausCoor[(2-1)*GausNum + 1-1] = 0.0;
    GausWeight[1-1] = 0.04166667;
    
    GausCoor[(1-1)*GausNum + 2-1] = 0.5;
    GausCoor[(2-1)*GausNum + 2-1] = 0.5;
    GausWeight[2-1] = 0.12500000;
    
    GausCoor[(1-1)*GausNum + 3-1] = 0.0;
    GausCoor[(2-1)*GausNum + 3-1] = 1.0;
    GausWeight[3-1] = 0.04166667;
    
    GausCoor[(1-1)*GausNum + 4-1] = 0.0;
    GausCoor[(2-1)*GausNum + 4-1] = 0.5;
    GausWeight[4-1] = 0.12500000;
    
    GausCoor[(1-1)*GausNum + 5-1] = 0.0;
    GausCoor[(2-1)*GausNum + 5-1] = 0.0;
    GausWeight[5-1] = 0.04166667;
    
    GausCoor[(1-1)*GausNum + 6-1] = 0.5;
    GausCoor[(2-1)*GausNum + 6-1] = 0.0;
    GausWeight[6-1] = 0.12500000;
    
    return;
}
void SetQuadrIntGaus4()  //Q
{
    GausNum = 4;
    
    GausCoor[(1-1)*GausNum + 1-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 1-1] =  5.773502692e-001;
    GausWeight[1-1] = 1.000000000e+000;
    
    GausCoor[(1-1)*GausNum + 2-1] =  5.773502692e-001;
    GausCoor[(2-1)*GausNum + 2-1] = -5.773502692e-001;
    GausWeight[2-1] = 1.000000000e+000;
    
    GausCoor[(1-1)*GausNum + 3-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 3-1] =  5.773502692e-001;
    GausWeight[3-1] = 1.000000000e+000;
    
    GausCoor[(1-1)*GausNum + 4-1] = -5.773502692e-001;
    GausCoor[(2-1)*GausNum + 4-1] = -5.773502692e-001;
    GausWeight[4-1] = 1.000000000e+000;
    
    return;
}
void SetQuadrIntGaus9()  //Q
{
    GausNum = 9;
    
    GausCoor[(1-1)*GausNum + 1-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum + 1-1] =  7.745966692e-001;
    GausWeight[1-1] = 3.086419753e-001;
    
    GausCoor[(1-1)*GausNum + 2-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum + 2-1] =  0.000000000e+000;
    GausWeight[2-1] = 4.938271605e-001;
    
    GausCoor[(1-1)*GausNum + 3-1] =  7.745966692e-001;
    GausCoor[(2-1)*GausNum + 3-1] = -7.745966692e-001;
    GausWeight[3-1] = 3.086419753e-001;
    
    GausCoor[(1-1)*GausNum + 4-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 4-1] =  7.745966692e-001;
    GausWeight[4-1] = 4.938271605e-001;
    
    GausCoor[(1-1)*GausNum + 5-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 5-1] =  0.000000000e+000;
    GausWeight[5-1] = 7.901234568e-001;
    
    GausCoor[(1-1)*GausNum + 6-1] =  0.000000000e+000;
    GausCoor[(2-1)*GausNum + 6-1] = -7.745966692e-001;
    GausWeight[6-1] = 4.938271605e-001;
    
    GausCoor[(1-1)*GausNum + 7-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 7-1] =  7.745966692e-001;
    GausWeight[7-1] = 3.086419753e-001;
    
    GausCoor[(1-1)*GausNum + 8-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 8-1] =  0.000000000e+000;
    GausWeight[8-1] = 4.938271605e-001;
    
    GausCoor[(1-1)*GausNum + 9-1] = -7.745966692e-001;
    GausCoor[(2-1)*GausNum + 9-1] = -7.745966692e-001;
    GausWeight[9-1] = 3.086419753e-001;
    
    return;
}
void SetQuadrIntNode4()  //Q
{
    GausNum = 4;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausCoor[(2-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 2-1] =  1.0;
    GausCoor[(2-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 3-1] =  1.0;
    GausCoor[(2-1)*GausNum + 3-1] =  1.0;
    GausWeight[3-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 4-1] = -1.0;
    GausCoor[(2-1)*GausNum + 4-1] =  1.0;
    GausWeight[4-1] = 1.0;
    
    return;
}
void SetQuadrIntNode8()  //Q
{
    GausNum = 8;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausCoor[(2-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 1.0/3.0;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.0;
    GausCoor[(2-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 2.0/3.0;
    
    GausCoor[(1-1)*GausNum + 3-1] =  1.0;
    GausCoor[(2-1)*GausNum + 3-1] = -1.0;
    GausWeight[3-1] = 1.0/3.0;
    
    GausCoor[(1-1)*GausNum + 4-1] =  1.0;
    GausCoor[(2-1)*GausNum + 4-1] =  0.0;
    GausWeight[4-1] = 2.0/3.0;
    
    GausCoor[(1-1)*GausNum + 5-1] =  1.0;
    GausCoor[(2-1)*GausNum + 5-1] =  1.0;
    GausWeight[5-1] = 1.0/3.0;
    
    GausCoor[(1-1)*GausNum + 6-1] =  0.0;
    GausCoor[(2-1)*GausNum + 6-1] =  1.0;
    GausWeight[6-1] = 2.0/3.0;
    
    GausCoor[(1-1)*GausNum + 7-1] = -1.0;
    GausCoor[(2-1)*GausNum + 7-1] =  1.0;
    GausWeight[7-1] = 1.0/3.0;
    
    GausCoor[(1-1)*GausNum + 8-1] = -1.0;
    GausCoor[(2-1)*GausNum + 8-1] =  0.0;
    GausWeight[8-1] = 2.0/3.0;
    
    return;
}
void SetQuadrIntNode9()  //Q
{
    GausNum = 9;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausCoor[(2-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.0;
    GausCoor[(2-1)*GausNum + 2-1] = -1.0;
    GausWeight[2-1] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 3-1] =  1.0;
    GausCoor[(2-1)*GausNum + 3-1] = -1.0;
    GausWeight[3-1] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 4-1] = -1.0;
    GausCoor[(2-1)*GausNum + 4-1] =  0.0;
    GausWeight[4-1] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 5-1] =  0.0;
    GausCoor[(2-1)*GausNum + 5-1] =  0.0;
    GausWeight[5-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 6-1] =  1.0;
    GausCoor[(2-1)*GausNum + 6-1] =  0.0;
    GausWeight[6-1] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 7-1] = -1.0;
    GausCoor[(2-1)*GausNum + 7-1] =  1.0;
    GausWeight[7-1] = 2.500000e-01;
    
    GausCoor[(1-1)*GausNum + 8-1] =  0.0;
    GausCoor[(2-1)*GausNum + 8-1] =  1.0;
    GausWeight[8-1] = 5.000000e-01;
    
    GausCoor[(1-1)*GausNum + 9-1] =  1.0;
    GausCoor[(2-1)*GausNum + 9-1] =  1.0;
    GausWeight[9-1] = 2.500000e-01;
    
    return;
}
void SetLineIntGaus2()   //L
{
    GausNum = 2;
    
    GausCoor[(1-1)*GausNum + 1-1] =  5.773502692e-001;
    GausWeight[1-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 2-1] = -5.773502692e-001;
    GausWeight[2-1] = 1.0;
    
    return;
}
void SetLineIntGaus3()   //L
{
    GausNum = 3;
    
    GausCoor[(1-1)*GausNum + 1-1] =  7.745966692e-001;
    GausWeight[1-1] = 5.555555556e-001;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.000000000e+000;
    GausWeight[2-1] = 8.888888889e-001;
    
    GausCoor[(1-1)*GausNum + 3-1] = -7.745966692e-001;
    GausWeight[3-1] = 5.555555556e-001;
    
    return;
}
void SetLineIntNode2()   //L
{
    GausNum = 2;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 2-1] =  1.0;
    GausWeight[2-1] = 1.0;
    
    return;
}
void SetLineIntNode3()   //L
{
    GausNum = 3;
    
    GausCoor[(1-1)*GausNum + 1-1] = -1.0;
    GausWeight[1-1] = 0.5;
    
    GausCoor[(1-1)*GausNum + 2-1] =  0.0;
    GausWeight[2-1] = 1.0;
    
    GausCoor[(1-1)*GausNum + 3-1] =  1.0;
    GausWeight[3-1] = 0.5;
    
    return;
}
void SetPointIntNode1()  //P
{
    GausNum = 1;
    
    GausCoor[(1-1)*GausNum + 1-1] = 0.0;
    GausWeight[1-1] = 1.0;
    
    return;
}




