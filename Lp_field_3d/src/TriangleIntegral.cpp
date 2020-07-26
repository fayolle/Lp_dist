#include <cassert>
#include <stdexcept>

#include "TriangleIntegral.h"
#include "GeometryUtils.h"


// Convert from int to enum and enum to int
template<>
NumSamples my_enum_convert<NumSamples>(int in) {
  switch(in) {
    case 1: return one_point;
    case 3: return three_points;
    case 7: return seven_points;
    case 24: return twenty_four_points;
    case 27: return twenty_seven_points;
    case 32: return thirty_two_points;
    default: throw std::logic_error(__FILE__ ": enum NumSamples out of range");
  }
}

/*
  int my_enum_convert(NumSamples in) {
  switch(in) {
    case one_point: return 1;
    case three_points: return 3;
    case seven_points: return 7;
    case twenty_four_points: return 24;
    case twenty_seven_points: return 27;
    case thirty_two_points: return 32;
  }
}
*/


const double TriangleIntegral::_samples_weights_1[] = {
  1.0/3, 1.0/3, 1.0
};

const double TriangleIntegral::_samples_weights_3[] = {
  1.0/6, 1.0/6, 1.0/3,
  1.0/6, 2.0/3, 1.0/3,
  2.0/3, 1.0/6, 1.0/3
};

const double TriangleIntegral::_samples_weights_4[] = {
  1.0/3, 1.0/3, -27.0/48,
  0.6, 0.2, 25.0/48,
  0.2, 0.6, 25.0/48,
  0.2, 0.2, 25.0/48
};

const double TriangleIntegral::_samples_weights_7[] = {
  1.0/3, 1.0/3, 0.225,
  0.0597158717, 0.4701420641, 0.1323941527,
  0.4701420641, 0.0597158717, 0.1323941527,
  0.4701420641, 0.4701420641, 0.1323941527,
  0.7974269853, 0.1012865073, 0.1259391805,
  0.1012865073, 0.7974269853, 0.1259391805,
  0.1012865073, 0.1012865073, 0.1259391805
};

const double TriangleIntegral::_samples_weights_24[] = {
  5.0550507373529086e-01, 2.0776116575484826e-01, 1.7344807725532943e-01,
  2.7542385024412980e-01, 4.8123289062464247e-01, 1.9053311454269983e-01,
  2.6481531651496770e-01, 2.7586334089315967e-01, 1.6882888511942015e-01,
  7.5329402776254240e-01, 1.0954959855585467e-01, 1.0546076281767805e-01,
  5.2433682558924433e-01, 3.6419744430339263e-01, 1.4815929467355968e-01,
  2.9530445535851102e-01, 6.4203365318662664e-01, 1.0983120878770872e-01,
  1.0614642990289996e-01, 7.6777680170023954e-01, 1.0507331820482332e-01,
  6.3491832379200652e-01, 3.6036266787907723e-02, 8.5924658784158670e-02,
  3.8729657913960353e-01, 8.4198522115543739e-02, 1.2537585060182724e-01,
  1.6929927488966462e-01, 1.0999439055630450e-01, 1.1594828119739846e-01,
  8.0491894656105567e-02, 5.7966325105486349e-01, 1.3237226895051976e-01,
  9.5379208487721689e-02, 3.3947290311800554e-01, 1.2348449173239080e-01,
  9.2899486985787905e-01, 4.7768381772022417e-02, 2.9216658446243379e-02,
  7.4726591728868819e-01, 2.2376358774275851e-01, 6.4605204046914597e-02,
  5.0365825075943971e-01, 4.8798437805397499e-01, 3.9118824435043810e-02,
  1.6134650499890957e-01, 8.3865349500109043e-01, 2.2133893564494179e-02,
  2.9553592846822851e-02, 9.3049846900263089e-01, 3.0406188052025412e-02,
  8.6854386943076545e-01, 3.8102570854643414e-03, 2.1333382551825181e-02,
  3.9366774470722010e-01, 0.0000000000000000e+00, 2.3800609628471206e-02,
  1.7690730625559031e-01, 1.0939142057119933e-02, 2.9693247293360987e-02,
  3.5319656252586096e-02, 3.9099745550423282e-02, 3.5311689185924387e-02,
  0.0000000000000000e+00, 7.7757518429429107e-01, 2.6798161571713618e-02,
  0.0000000000000000e+00, 4.6374383867430541e-01, 3.0312523835131357e-02,
  3.0573404093099332e-02, 1.9305903224251936e-01, 6.2829404721337689e-02
};

const double TriangleIntegral::_samples_weights_27[] = {
  4.6494564773693992e-01, 2.9133859436942361e-01, 1.3648275991498204e-01,
  3.2081957909482994e-01, 5.3634228112084714e-01, 1.2438630022250971e-01,
  5.1353143433447235e-01, 1.2454405910544103e-01, 1.1329177024539897e-01,
  2.8790310224819649e-01, 2.2789955884347501e-01, 1.3228489176992250e-01,
  2.6677168071577745e-01, 4.1132499178904658e-01, 1.1722353681481934e-01,
  1.1698976413323442e-01, 3.1909737814681871e-01, 1.0998202543484477e-01,
  8.1626233715968810e-01, 2.7719522918618567e-02, 4.7284119131529377e-02,
  5.6938486195327997e-01, 3.4992914334288650e-01, 1.0994399601768742e-01,
  3.7272769861629096e-01, 5.9895439629934211e-01, 6.5193746289815974e-02,
  2.6807150626772580e-02, 8.1562969693268217e-01, 4.6224760707242137e-02,
  7.0099267949645228e-01, 1.4118119730952799e-01, 1.0412107067624195e-01,
  3.2719878157552895e-01, 8.1721404855381763e-02, 8.5195409796230526e-02,
  1.3667083534390506e-01, 1.3035453031942690e-01, 9.1076518240300441e-02,
  1.3828000204292318e-01, 7.1027868107761583e-01, 9.8381989816749074e-02,
  2.2592651051306589e-02, 3.8913981113319357e-01, 5.3445574349465230e-02,
  9.3614893514675623e-01, 3.2899822292186298e-02, 2.6211869704176473e-02,
  8.0454974747615537e-01, 1.6429286715713465e-01, 5.5191800300359820e-02,
  6.1948431533135195e-01, 3.7802163891336921e-01, 2.2550142431420638e-02,
  1.6655614492060572e-01, 8.0364834053903877e-01, 5.3513272326506316e-02,
  3.3268560622678411e-02, 9.3551434285897095e-01, 2.6748618572925459e-02,
  6.1924873232110123e-01, 2.6297199713764152e-02, 5.8869116212867049e-02,
  3.9659731669586495e-01, 1.4354532010930898e-02, 3.6717768780272685e-02,
  1.6892970982290229e-01, 2.2120535196161750e-02, 4.2755616195827365e-02,
  3.2916403878999745e-02, 3.4222771841359190e-02, 2.9096217361124159e-02,
  2.5660186833052434e-02, 6.1758873171277151e-01, 5.7443554735054178e-02,
  1.2417148586801485e-01, 5.3141960154079959e-01, 1.0824295295050959e-01,
  2.5252704638304480e-02, 1.7400571673032256e-01, 4.8140601001216463e-02
};

const double TriangleIntegral::_samples_weights_32[] = {
  3.7986021093401956e-01, 2.1078525939140391e-01, 1.1887566790227083e-01,
  3.0141709320909305e-01, 4.0978657777002531e-01, 1.5044412520664885e-01,
  5.5802528953120256e-01, 2.1377743253005960e-01, 1.2632909284531338e-01,
  1.2512299505810387e-01, 6.1938125736255578e-01, 1.0192984975357525e-01,
  2.1117939909804934e-01, 2.4498296509349016e-01, 9.4999150650614317e-02,
  8.5431474947580432e-01, 7.1871496101589105e-02, 4.4981492398316447e-02,
  7.1788185898052326e-01, 2.0376848107772977e-01, 7.9147211585943858e-02,
  4.6631787462323071e-01, 4.0896380449124475e-01, 1.1997941465421234e-01,
  2.5015500335339214e-01, 6.2768261568031403e-01, 1.0670416609764186e-01,
  7.9955384841381316e-02, 8.2600331401756000e-01, 6.1058344824144795e-02,
  7.1008125956836521e-01, 6.4413220382260550e-02, 8.2563774790925248e-02,
  4.9732063377796598e-01, 7.0566724344036824e-02, 9.6297610073814668e-02,
  2.6077068256562896e-01, 9.5428585810584610e-02, 9.1875684331583440e-02,
  8.9602705800587434e-02, 1.1638649906727733e-01, 6.1150555208077911e-02,
  2.3088148766115757e-02, 7.4918973979067949e-01, 4.3370170834023010e-02,
  1.2953296900433620e-01, 4.2260565743346001e-01, 1.0829374522633514e-01,
  9.3448087604440955e-02, 2.4345813394879973e-01, 5.5887468639759713e-02,
  9.5526919357006035e-01, 2.3551733249578710e-02, 1.3351800054734712e-02,
  8.4593539837314391e-01, 1.5406460162685609e-01, 1.5428984747249670e-02,
  6.1600929617267497e-01, 3.6118159118967208e-01, 5.0198346855370224e-02,
  3.9316510319604808e-01, 5.8168921474014745e-01, 5.6291117210426664e-02,
  1.8920633061715936e-01, 7.8860171922313160e-01, 4.1240008239364231e-02,
  4.3010560106405471e-02, 9.4547507322097091e-01, 1.4239502872161450e-02,
  8.5815888421533082e-01, 0.0000000000000000e+00, 1.3691069308687381e-02,
  6.2731531923241179e-01, 0.0000000000000000e+00, 1.9309417484872689e-02,
  3.6384660446077510e-01, 1.4566514788346974e-02, 3.7090960843213061e-02,
  1.5557066896897953e-01, 2.1152223383121949e-02, 3.6967371622461546e-02,
  2.9754117496841759e-02, 2.7110971356255786e-02, 2.1018653471205032e-02,
  0.0000000000000000e+00, 9.2734897448394982e-01, 9.7760996293200769e-03,
  2.5716283623693881e-02, 5.4444667627192522e-01, 5.6339308919459923e-02,
  2.4506286636990005e-02, 3.3212908394764507e-01, 4.9808146403015403e-02,
  9.2296909059649268e-03, 1.4604496167217568e-01, 2.1361687315256585e-02
};

TriangleIntegral::TriangleIntegral(NumSamples num_samples, const Point3& p1,
                                   const Point3& p2, const Point3& p3) {

  _num_samples = num_samples;
  const double* samples_weights = NULL;
  
  if (num_samples == 1) {
    _weights = new double[1];
    _samples = new Point3[1];

    samples_weights = _samples_weights_1;
    
  } else if (num_samples == 3) {
    _weights = new double[3];
    _samples = new Point3[3];

    samples_weights = _samples_weights_3;
    
  } else if (num_samples == 4) {
    _weights = new double[4];
    _samples = new Point3[4];

    samples_weights = _samples_weights_4;

  } else if (num_samples == 7) {
    _weights = new double[7];
    _samples = new Point3[7];

    samples_weights = _samples_weights_7;
    
  } else if (num_samples == 24) {
    _weights = new double[24];
    _samples = new Point3[24];

    samples_weights = _samples_weights_24;
    
  } else if (num_samples == 27) {
    _weights = new double[27];
    _samples = new Point3[27];

    samples_weights = _samples_weights_27;

  } else if (num_samples == 32) {
    _weights = new double[32];
    _samples = new Point3[32];

    samples_weights = _samples_weights_32;
    
  }

  //Point3 q1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
  Vector3 q1 = p2 - p1;
  //Point3 q2(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
  Vector3 q2 = p3 - p1;

  _area = norm( cross(q1, q2) ) / 2.0;
  
  for (unsigned int i = 0; i < _num_samples; ++i) {
    Point3 sample = p1 + q1 * samples_weights[3*i + 0] +
        q2 * samples_weights[3*i + 1];
    _samples[i] = sample;
    _weights[i] = samples_weights[3 * i + 2] * _area;
  }

  double weightsum = 0.0;
  for (unsigned int i = 0; i < _num_samples; ++i) {
    weightsum += samples_weights[3 * i + 2];
  }

  for (unsigned int i = 0; i < _num_samples; ++i) {
    _weights[i] /= weightsum;
  }
}

TriangleIntegral::~TriangleIntegral() {
  delete[] _weights;
  delete[] _samples;
}

Point3 TriangleIntegral::get_sample(int i) {
  return _samples[i];
}

double TriangleIntegral::get_weight(int i) {
  return _weights[i];
}

double TriangleIntegral::get_area() {
  return _area;
}

unsigned int TriangleIntegral::get_num_samples() {
  return _num_samples;
}

double TriangleIntegral::integrate(double fi[]) {
  double integral = 0;
  for (unsigned int i = 0; i < _num_samples; ++i) {
    integral += fi[i] * _weights[i];
  }
  return integral;
}

double TriangleIntegral::integrate(const std::vector<double>& fi) {
  double integral = 0.0;
  assert(_num_samples == fi.size());
  for (unsigned int i = 0; i < _num_samples; ++i) {
    integral += fi[i] * _weights[i];
  }
  return integral;
}

// Test for the class TriangleIntegral
#ifdef TEST
#include <iostream>
using std::cout;
using std::endl;

int main (void) {
  // Domain D is a square: [0,0] x [1,1]
  Point3 p1(0.0, 0.0, 0.0);
  Point3 p2(1.0, 0.0, 0.0);
  Point3 p3(0.0, 1.0, 0.0);
  Point3 p4(1.0, 1.0, 0.0);

  // Compute Int[x, D]
  NumSamples num_samples = three_points;
  TriangleIntegral t1(three_points, p1, p2, p3);
  cout << "Area = " << t1.get_area() << endl;
  
  double fi[3];
  for (int i = 0; i < t1.get_num_samples(); i++) {
    Point3 s = t1.get_sample(i);
    fi[i] = s.x;
    cout << "fi[" << i << "] = " << s.x << endl;
  }
  double int1 = t1.integrate(fi);
  cout << "int1 = " << int1 << endl;

  TriangleIntegral t2(num_samples, p2, p4, p3);
  cout << "Area = " << t2.get_area() << endl;
  
  for (int i = 0; i < t2.get_num_samples(); ++i) {
    Point3 s = t2.get_sample(i);
    fi[i] = s.x;
    cout << "fi[" << i << "] = " << s.x << endl;
  }
  double int2 = t2.integrate(fi);
  cout << "int2 = " << int2 << endl;

  double result = int1 + int2;
  cout << "Int[x, D] = " << result << endl;


  // Compute Int[1, D]
  for (int i = 0; i < t1.get_num_samples(); i++) {
    fi[i] = 1.0;
    cout << "fi[" << i << "] = 1.0" << endl;
  }
  int1 = t1.integrate(fi);
  cout << "int1 = " << int1 << endl;

  for (int i = 0; i < t2.get_num_samples(); ++i) {
    fi[i] = 1.0;
    cout << "fi[" << i << "] = 1.0" << endl;
  }
  int2 = t2.integrate(fi);
  cout << "int2 = " << int2 << endl;

  result = int1 + int2;
  cout << "Int[1, D] = " << result << endl;


  // Compute Int[x, D2] where D2 is [0,0]x[1/2,1/2]
  Point3 q1(0.0, 0.0, 0.0);
  Point3 q2(0.5, 0.0, 0.0);
  Point3 q3(0.0, 0.5, 0.0);
  Point3 q4(0.5, 0.5, 0.0);

  // Compute Int[x, D]
  num_samples = three_points;
  TriangleIntegral tt1(num_samples, q1, q2, q3);
  cout << "Area = " << tt1.get_area() << endl;
  
  for (int i = 0; i < tt1.get_num_samples(); i++) {
    Point3 s = tt1.get_sample(i);
    fi[i] = s.x;
    cout << "fi[" << i << "] = " << s.x << endl;
  }
  int1 = tt1.integrate(fi);
  cout << "int1 = " << int1 << endl;

  TriangleIntegral tt2(num_samples, q2, q4, q3);
  cout << "Area = " << tt2.get_area() << endl;
  
  for (int i = 0; i < tt2.get_num_samples(); ++i) {
    Point3 s = tt2.get_sample(i);
    fi[i] = s.x;
    cout << "fi[" << i << "] = " << s.x << endl;
  }
  int2 = tt2.integrate(fi);
  cout << "int2 = " << int2 << endl;

  result = int1 + int2;
  cout << "Int[x, D2] = " << result << endl;
  
  return 0;
}

#endif


// For speed: get rid of the class for the special cases of 27 and 32 samples
//
const double samples_weights_27[] = {
4.6494564773693992e-01, 2.9133859436942361e-01, 1.3648275991498204e-01,
3.2081957909482994e-01, 5.3634228112084714e-01, 1.2438630022250971e-01,
5.1353143433447235e-01, 1.2454405910544103e-01, 1.1329177024539897e-01,
2.8790310224819649e-01, 2.2789955884347501e-01, 1.3228489176992250e-01,
2.6677168071577745e-01, 4.1132499178904658e-01, 1.1722353681481934e-01,
1.1698976413323442e-01, 3.1909737814681871e-01, 1.0998202543484477e-01,
8.1626233715968810e-01, 2.7719522918618567e-02, 4.7284119131529377e-02,
5.6938486195327997e-01, 3.4992914334288650e-01, 1.0994399601768742e-01,
3.7272769861629096e-01, 5.9895439629934211e-01, 6.5193746289815974e-02,
2.6807150626772580e-02, 8.1562969693268217e-01, 4.6224760707242137e-02,
7.0099267949645228e-01, 1.4118119730952799e-01, 1.0412107067624195e-01,
3.2719878157552895e-01, 8.1721404855381763e-02, 8.5195409796230526e-02,
1.3667083534390506e-01, 1.3035453031942690e-01, 9.1076518240300441e-02,
1.3828000204292318e-01, 7.1027868107761583e-01, 9.8381989816749074e-02,
2.2592651051306589e-02, 3.8913981113319357e-01, 5.3445574349465230e-02,
9.3614893514675623e-01, 3.2899822292186298e-02, 2.6211869704176473e-02,
8.0454974747615537e-01, 1.6429286715713465e-01, 5.5191800300359820e-02,
6.1948431533135195e-01, 3.7802163891336921e-01, 2.2550142431420638e-02,
1.6655614492060572e-01, 8.0364834053903877e-01, 5.3513272326506316e-02,
3.3268560622678411e-02, 9.3551434285897095e-01, 2.6748618572925459e-02,
6.1924873232110123e-01, 2.6297199713764152e-02, 5.8869116212867049e-02,
3.9659731669586495e-01, 1.4354532010930898e-02, 3.6717768780272685e-02,
1.6892970982290229e-01, 2.2120535196161750e-02, 4.2755616195827365e-02,
3.2916403878999745e-02, 3.4222771841359190e-02, 2.9096217361124159e-02,
2.5660186833052434e-02, 6.1758873171277151e-01, 5.7443554735054178e-02,
1.2417148586801485e-01, 5.3141960154079959e-01, 1.0824295295050959e-01,
2.5252704638304480e-02, 1.7400571673032256e-01, 4.8140601001216463e-02
};

const double samples_weights_32[] = {
3.7986021093401956e-01, 2.1078525939140391e-01, 1.1887566790227083e-01,
3.0141709320909305e-01, 4.0978657777002531e-01, 1.5044412520664885e-01,
5.5802528953120256e-01, 2.1377743253005960e-01, 1.2632909284531338e-01,
1.2512299505810387e-01, 6.1938125736255578e-01, 1.0192984975357525e-01,
2.1117939909804934e-01, 2.4498296509349016e-01, 9.4999150650614317e-02,
8.5431474947580432e-01, 7.1871496101589105e-02, 4.4981492398316447e-02,
7.1788185898052326e-01, 2.0376848107772977e-01, 7.9147211585943858e-02,
4.6631787462323071e-01, 4.0896380449124475e-01, 1.1997941465421234e-01,
2.5015500335339214e-01, 6.2768261568031403e-01, 1.0670416609764186e-01,
7.9955384841381316e-02, 8.2600331401756000e-01, 6.1058344824144795e-02,
7.1008125956836521e-01, 6.4413220382260550e-02, 8.2563774790925248e-02,
4.9732063377796598e-01, 7.0566724344036824e-02, 9.6297610073814668e-02,
2.6077068256562896e-01, 9.5428585810584610e-02, 9.1875684331583440e-02,
8.9602705800587434e-02, 1.1638649906727733e-01, 6.1150555208077911e-02,
2.3088148766115757e-02, 7.4918973979067949e-01, 4.3370170834023010e-02,
1.2953296900433620e-01, 4.2260565743346001e-01, 1.0829374522633514e-01,
9.3448087604440955e-02, 2.4345813394879973e-01, 5.5887468639759713e-02,
9.5526919357006035e-01, 2.3551733249578710e-02, 1.3351800054734712e-02,
8.4593539837314391e-01, 1.5406460162685609e-01, 1.5428984747249670e-02,
6.1600929617267497e-01, 3.6118159118967208e-01, 5.0198346855370224e-02,
3.9316510319604808e-01, 5.8168921474014745e-01, 5.6291117210426664e-02,
1.8920633061715936e-01, 7.8860171922313160e-01, 4.1240008239364231e-02,
4.3010560106405471e-02, 9.4547507322097091e-01, 1.4239502872161450e-02,
8.5815888421533082e-01, 0.0000000000000000e+00, 1.3691069308687381e-02,
6.2731531923241179e-01, 0.0000000000000000e+00, 1.9309417484872689e-02,
3.6384660446077510e-01, 1.4566514788346974e-02, 3.7090960843213061e-02,
1.5557066896897953e-01, 2.1152223383121949e-02, 3.6967371622461546e-02,
2.9754117496841759e-02, 2.7110971356255786e-02, 2.1018653471205032e-02,
0.0000000000000000e+00, 9.2734897448394982e-01, 9.7760996293200769e-03,
2.5716283623693881e-02, 5.4444667627192522e-01, 5.6339308919459923e-02,
2.4506286636990005e-02, 3.3212908394764507e-01, 4.9808146403015403e-02,
9.2296909059649268e-03, 1.4604496167217568e-01, 2.1361687315256585e-02
};

inline double computeArea(const Point3& p1, const Point3& p2, 
						  const Point3& p3) {
return norm(cross(p2-p1,p3-p1)) / 2.0;
}

void initArrays27(const Point3& p1, const Point3& p2, const Point3& p3,
				Point3 samples[], double weights[]) {
	double area = computeArea(p1, p2, p3);
	for (unsigned int i = 0; i < 27; ++i) {
		samples[i] = p1 + (p2 - p1)*samples_weights_27[3*i + 0] +
        (p3 - p1)*samples_weights_27[3*i + 1];
		weights[i] = samples_weights_27[3*i + 2] * area;
	}
	double weightsum = 0.0;
	for (unsigned int i = 0; i < 27; ++i) {
		weightsum += samples_weights_27[3*i + 2];
	}
	for (unsigned int i = 0; i < 27; ++i) {
		weights[i] /= weightsum;
	}
}

void initArrays32(const Point3& p1, const Point3& p2, const Point3& p3,
				Point3 samples[], double weights[]) {
	double area = computeArea(p1, p2, p3);
	for (unsigned int i = 0; i < 32; ++i) {
		samples[i] = p1 + (p2 - p1)*samples_weights_32[3*i + 0] +
        (p3 - p1)*samples_weights_32[3*i + 1];
		weights[i] = samples_weights_32[3*i + 2] * area;
	}
	double weightsum = 0.0;
	for (unsigned int i = 0; i < 32; ++i) {
		weightsum += samples_weights_32[3*i + 2];
	}
	for (unsigned int i = 0; i < 32; ++i) {
		weights[i] /= weightsum;
	}
}

