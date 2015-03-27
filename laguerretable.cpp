#include "./classes.h"

void laguerretable(int n, rl* glx, rl* glw){

  //Rather than recompute gauss legendre points and weights every time they
  //are needed, it will be faster to precompute them for small n.
  // Here, values for n<=50 were calculated using scipy

static rl x1[1]={1.000000000000000e+00};
static rl w1[1]={1.000000000000000e+00};

static rl x2[2]={5.857864376269051e-01, 3.414213562373095e+00};
static rl w2[2]={8.535533905932727e-01, 1.464466094067263e-01};

static rl x3[3]={4.157745567834792e-01, 2.294280360279040e+00, 6.289945082937479e+00};
static rl w3[3]={7.110930099291716e-01, 2.785177335692420e-01, 1.038925650158615e-02};

static rl x4[4]={3.225476896193926e-01, 1.745761101158347e+00, 4.536620296921130e+00, 9.395070912301133e+00};
static rl w4[4]={6.031541043416266e-01, 3.574186924377983e-01, 3.888790851500527e-02, 5.392947055613273e-04};

static rl x5[5]={2.635603197181406e-01, 1.413403059106515e+00, 3.596425771040719e+00, 7.085810005858837e+00, 1.264080084427578e+01};
static rl w5[5]={5.217556105828127e-01, 3.986668110831817e-01, 7.594244968170824e-02, 3.611758679922049e-03, 2.336997238577623e-05};

static rl x6[6]={2.228466041792603e-01, 1.188932101672623e+00, 2.992736326059316e+00, 5.775143569104514e+00, 9.837467418382598e+00, 1.598287398060170e+01};
static rl w6[6]={4.589646739499807e-01, 4.170008307721178e-01, 1.133733820740440e-01, 1.039919745314900e-02, 2.610172028149294e-04, 8.985479064296205e-07};

static rl x7[7]={1.930436765603619e-01, 1.026664895339191e+00, 2.567876744950746e+00, 4.900353084526482e+00, 8.182153444562864e+00, 1.273418029179782e+01, 1.939572786226254e+01};
static rl w7[7]={4.093189517012897e-01, 4.218312778617227e-01, 1.471263486575051e-01, 2.063351446871710e-02, 1.074010143280738e-03, 1.586546434856411e-05, 3.170315478995584e-08};

static rl x8[8]={1.702796323051006e-01, 9.037017767993795e-01, 2.251086629866132e+00, 4.266700170287656e+00, 7.045905402393460e+00, 1.075851601018099e+01, 1.574067864127800e+01, 2.286313173688926e+01};
static rl w8[8]={3.691885893416525e-01, 4.187867808143441e-01, 1.757949866371706e-01, 3.334349226121599e-02, 2.794536235225710e-03, 9.076508773358259e-05, 8.485746716272535e-07, 1.048001174871513e-09};

static rl x9[9]={1.523222277318079e-01, 8.072200227422557e-01, 2.005135155619346e+00, 3.783473973331232e+00, 6.204956777876609e+00, 9.372985251687572e+00, 1.346623691109209e+01, 1.883359778899171e+01, 2.637407189092737e+01};
static rl w9[9]={3.361264217979638e-01, 4.112139804239895e-01, 1.992875253708887e-01, 4.746056276565175e-02, 5.599626610794632e-03, 3.052497670932133e-04, 6.592123026075409e-06, 4.110769330349489e-08, 3.290874030350722e-11};

static rl x10[10]={1.377934705404925e-01, 7.294545495031682e-01, 1.808342901740315e+00, 3.401433697854900e+00, 5.552496140063803e+00, 8.330152746764485e+00, 1.184378583790007e+01, 1.627925783137810e+01, 2.199658581198078e+01, 2.992069701227389e+01};
static rl w10[10]={3.084411157650263e-01, 4.011199291552968e-01, 2.180682876118141e-01, 6.208745609867743e-02, 9.501516975181139e-03, 7.530083885875612e-04, 2.825923349599526e-05, 4.249313984962677e-07, 1.839564823979602e-09, 9.911827219609021e-13};

static rl x11[11]={1.257964421879685e-01, 6.654182558392288e-01, 1.647150545872170e+00, 3.091138143035250e+00, 5.029284401579837e+00, 7.509887863806623e+00, 1.060595099954698e+01, 1.443161375806420e+01, 1.917885740321466e+01, 2.521770933967756e+01, 3.349719284717554e+01};
static rl w11[11]={2.849332128941453e-01, 3.897208895278430e-01, 2.327818318489908e-01, 7.656445354619940e-02, 1.439328276735048e-02, 1.518880846484847e-03, 8.513122435471719e-05, 2.292403879574458e-06, 2.486353702767852e-08, 7.712626933691337e-11, 2.883775868323619e-14};

static rl x12[12]={1.157221173580204e-01, 6.117574845151292e-01, 1.512610269776415e+00, 2.833751337743510e+00, 4.599227639418352e+00, 6.844525453115179e+00, 9.621316842456860e+00, 1.300605499330636e+01, 1.711685518746227e+01, 2.215109037939699e+01, 2.848796725098403e+01, 3.709912104446693e+01};
static rl w12[12]={2.647313710554613e-01, 3.777592758731565e-01, 2.440820113198909e-01, 9.044922221167902e-02, 2.010238115463374e-02, 2.663973541865301e-03, 2.032315926630033e-04, 8.365055856819621e-06, 1.668493876540864e-07, 1.342391030515028e-09, 3.061601635034946e-12, 8.148077467426205e-16};

static rl x13[13]={1.071423884722528e-01, 5.661318990404003e-01, 1.398564336451021e+00, 2.616597108406414e+00, 4.238845929017033e+00, 6.292256271140085e+00, 8.815001941186987e+00, 1.186140358881125e+01, 1.551076203770373e+01, 1.988463566388024e+01, 2.518526386467773e+01, 3.180038630194725e+01, 4.072300866926558e+01};
static rl w13[13]={2.471887084299251e-01, 3.656888229005547e-01, 2.525624200576507e-01, 1.034707580241813e-01, 2.643275441556180e-02, 4.220396040254564e-03, 4.118817704727242e-04, 2.351547398155296e-05, 7.317311620249345e-07, 1.108841625703970e-08, 6.770826692206133e-11, 1.159979959905098e-13, 2.245093203892754e-17};

static rl x14[14]={9.974750703259661e-02, 5.268576488519030e-01, 1.300629121251496e+00, 2.430801078730845e+00, 3.932102822293218e+00, 5.825536218301706e+00, 8.140240141565155e+00, 1.091649950736602e+01, 1.421080501116129e+01, 1.810489222021808e+01, 2.272338162826963e+01, 2.827298172324823e+01, 3.514944366059243e+01, 4.436608171111737e+01};
static rl w14[14]={2.318155771449443e-01, 3.537846915975388e-01, 2.587346102454318e-01, 1.154828935569227e-01, 3.319209215933770e-02, 6.192869437006677e-03, 7.398903778673614e-04, 5.490719466841715e-05, 2.409585764085360e-06, 5.801543981676702e-08, 6.819314692484947e-10, 3.221207751894758e-12, 4.221352440516562e-15, 6.052375022289402e-19};

static rl x15[15]={9.330781201728239e-02, 4.926917403018868e-01, 1.215595412070952e+00, 2.269949526203738e+00, 3.667622721751438e+00, 5.425336627413562e+00, 7.565916226613071e+00, 1.012022856801910e+01, 1.313028248217574e+01, 1.665440770832996e+01, 2.077647889944875e+01, 2.562389422672877e+01, 3.140751916975396e+01, 3.853068330648601e+01, 4.802608557268580e+01};
static rl w15[15]={2.182348859400322e-01, 3.422101779228286e-01, 2.630275779416589e-01, 1.264258181059394e-01, 4.020686492100059e-02, 8.563877803611453e-03, 1.212436147214238e-03, 1.116743923442548e-04, 6.459926762022725e-06, 2.226316907096268e-07, 4.227430384979485e-09, 3.921897267041128e-11, 1.456515264073096e-13, 1.483027051113300e-16, 1.600594906211127e-20};

static rl x16[16]={8.764941047892821e-02, 4.626963289150810e-01, 1.141057774831228e+00, 2.129283645098382e+00, 3.437086633893208e+00, 5.078018614549773e+00, 7.070338535048228e+00, 9.438314336391940e+00, 1.221422336886615e+01, 1.544152736878162e+01, 1.918015685675313e+01, 2.351590569399190e+01, 2.857872974288216e+01, 3.458339870228665e+01, 4.194045264768834e+01, 5.170116033954332e+01};
static rl w16[16]={2.061517149577735e-01, 3.310578549509001e-01, 2.657957776442049e-01, 1.362969342963740e-01, 4.732892869412422e-02, 1.129990008033912e-02, 1.849070943526365e-03, 2.042719153082775e-04, 1.484458687398155e-05, 6.828319330871196e-07, 1.881024841079679e-08, 2.862350242973900e-10, 2.127079033224054e-12, 6.297967002517717e-15, 5.050473700035486e-18, 4.161462370372846e-22};

static rl x17[17]={8.263821470894844e-02, 4.361503235587094e-01, 1.075176577511430e+00, 2.005193531649235e+00, 3.234256124047448e+00, 4.773513513700192e+00, 6.637829205364956e+00, 8.846685511169802e+00, 1.142552931937333e+01, 1.440782303748131e+01, 1.783828473070109e+01, 2.177826825772226e+01, 2.631531781124881e+01, 3.158177168045672e+01, 3.779609383747707e+01, 4.537571653398896e+01, 5.538975178983960e+01};
static rl w17[17]={1.953322052517190e-01, 3.203753572745567e-01, 2.673297263571680e-01, 1.451298543587538e-01, 5.443694324533601e-02, 1.435729776606245e-02, 2.662824735572731e-03, 3.436797271562963e-04, 3.027551783783043e-05, 1.768515053231715e-06, 6.576272886811079e-08, 1.469730932159567e-09, 1.816910362555415e-11, 1.095401388928706e-13, 2.617373882223454e-16, 1.672935693146166e-19, 1.065626316274051e-23};

static rl x18[18]={7.816916666970752e-02, 4.124900852591297e-01, 1.016520179623539e+00, 1.894888509969758e+00, 3.054353113202661e+00, 4.504205538889903e+00, 6.256725073949109e+00, 8.327825156605632e+00, 1.073799004775760e+01, 1.351365620755507e+01, 1.668930628193011e+01, 2.031076762626776e+01, 2.444068135928370e+01, 2.916820866257960e+01, 3.462792706566020e+01, 4.104181677280877e+01, 4.883392271608655e+01, 5.909054643590132e+01};
static rl w18[18]={1.855886031467488e-01, 3.101817663702299e-01, 2.678665671485377e-01, 1.529797474680835e-01, 6.143491786096070e-02, 1.768721308077147e-02, 3.660179767759964e-03, 5.406227870077285e-04, 5.616965051214407e-05, 4.015307883701343e-06, 1.914669856675652e-07, 5.836095268631428e-09, 1.071711266955398e-10, 1.089098713888859e-12, 5.386664748378172e-15, 1.049865978035697e-17, 5.405398451630925e-21, 2.691653269200917e-25};

static rl x19[19]={7.415878375719648e-02, 3.912686133199992e-01, 9.639573439979503e-01, 1.796175582068335e+00, 2.893651381873783e+00, 4.264215539627772e+00, 5.918141561644045e+00, 7.868618915334743e+00, 1.013242371681527e+01, 1.273088146384240e+01, 1.569127833983590e+01, 1.904899320982352e+01, 2.285084976082950e+01, 2.716066932741144e+01, 3.206912225186224e+01, 3.771290580121966e+01, 4.431736279583151e+01, 5.231290245740438e+01, 6.280242315350043e+01};
static rl w19[19]={1.767684749167235e-01, 3.004781436071457e-01, 2.675995470382572e-01, 1.599133721355544e-01, 6.824937997615102e-02, 2.123930760654337e-02, 4.841627351148501e-03, 8.049127473813304e-04, 9.652472093153445e-05, 8.207305258050858e-06, 4.830566724730696e-07, 1.904991361123398e-08, 4.816684630927940e-10, 7.348258839551273e-12, 6.202275387572637e-14, 2.541430843015391e-16, 4.078861296825672e-19, 1.707750187593836e-22, 6.715064649907943e-27};

static rl x20[20]={7.053988969198671e-02, 3.721268180016177e-01, 9.165821024832752e-01, 1.707306531028346e+00, 2.749199255309432e+00, 4.048925313850892e+00, 5.615174970861620e+00, 7.459017453671057e+00, 9.594392869581092e+00, 1.203880254696432e+01, 1.481429344263073e+01, 1.794889552051938e+01, 2.147878824028503e+01, 2.545170279318692e+01, 2.993255463170064e+01, 3.501343424047900e+01, 4.083305705672859e+01, 4.761999404734648e+01, 5.581079575006389e+01, 6.652441652561583e+01};
static rl w20[20]={1.687468018513514e-01, 2.912543620058509e-01, 2.666861028669822e-01, 1.660024532694981e-01, 7.482606466879266e-02, 2.496441730928206e-02, 6.202550844572084e-03, 1.144962386476953e-03, 1.557417730278146e-04, 1.540144086522482e-05, 1.086486366518023e-06, 5.330120909556676e-08, 1.757981179050528e-09, 3.725502402512247e-11, 4.767529251578026e-13, 3.372844243362437e-15, 1.155014339500373e-17, 1.539522140582371e-20, 5.286442725569206e-24, 1.656456612498941e-28};

static rl x21[21]={6.725781792316127e-02, 3.547728953235072e-01, 8.736601667786409e-01, 1.626869941929206e+00, 2.618626410545553e+00, 3.854652138109756e+00, 5.342369280622432e+00, 7.091168813219666e+00, 9.112778854269701e+00, 1.142177176237839e+01, 1.403627069787377e+01, 1.697895269278381e+01, 2.027850941499374e+01, 2.397184558715149e+01, 2.810752860094429e+01, 3.275149741056082e+01, 3.799718781926101e+01, 4.398524575714224e+01, 5.094735118993994e+01, 5.932599412023001e+01, 7.025568862801899e+01};
static rl w21[21]={1.614201001816523e-01, 2.824935645831712e-01, 2.652545755139150e-01, 1.713193963900418e-01, 8.112662589234178e-02, 2.881630176777059e-02, 7.734239426118855e-03, 1.567481390421274e-03, 2.384226587634183e-04, 2.693007343968514e-05, 2.224790793644991e-06, 1.317334497432410e-07, 5.443952290646474e-09, 1.516366205198029e-10, 2.717883118058966e-12, 2.942519354071984e-14, 1.759074035797611e-16, 5.073802227171921e-19, 5.659720870324156e-22, 1.606260924799283e-25, 4.044133927115794e-30};

static rl x22[22]={6.426762874480907e-02, 3.389672548149100e-01, 8.345899854491788e-01, 1.553713386750022e+00, 2.500006236736168e+00, 3.678420344463651e+00, 5.095349068968734e+00, 6.758835515813534e+00, 8.678853313678427e+00, 1.086768689354057e+01, 1.334045105149525e+01, 1.611581133876632e+01, 1.921700368490135e+01, 2.267331660143459e+01, 2.652231958763711e+01, 3.081335794450654e+01, 3.561333520319614e+01, 4.101696859847163e+01, 4.716675784915353e+01, 5.429738525657886e+01, 6.285709624629904e+01, 7.399550700859940e+01};
static rl w22[22]={1.547019876399348e-01, 2.741750823285843e-01, 2.634097015356312e-01, 1.759346092630954e-01, 8.712562217221634e-02, 3.275268789653316e-02, 9.424900256998582e-03, 2.077335752654816e-03, 3.491673714923000e-04, 4.438219314066161e-05, 4.214229901104530e-06, 2.940668246130056e-07, 1.476193826389758e-08, 5.186423025724909e-10, 1.230526251554477e-11, 1.880465509631537e-13, 1.735564139126912e-15, 8.831402308327667e-18, 2.160718646818175e-20, 2.031119967583367e-23, 4.797903520990772e-27, 9.780246506684377e-32};

static rl x23[23]={6.153203775752014e-02, 3.245113291545534e-01, 7.988739686929823e-01, 1.486886334222205e+00, 2.391755829792771e+00, 3.517797442503539e+00, 4.870560368733455e+00, 6.456991306767748e+00, 8.285651892597997e+00, 1.036700904344577e+01, 1.271382502816259e+01, 1.534168757288852e+01, 1.826974239389599e+01, 2.152172785081421e+01, 2.512747713011548e+01, 2.912517510875098e+01, 3.356489656311643e+01, 3.851445967718071e+01, 4.406980781487124e+01, 5.037522660304375e+01, 5.766830435043416e+01, 6.640287306834519e+01, 7.774322728471226e+01};
static rl w23[23]={1.485197990820924e-01, 2.662763422219940e-01, 2.612370674653731e-01, 1.799149140382408e-01, 9.280788216398320e-02, 3.673579259475661e-02, 1.126058922254327e-02, 2.677540143126828e-03, 4.923896728665736e-04, 6.956649894243718e-05, 7.476361523134264e-06, 6.030632412462762e-07, 3.588339690021517e-08, 1.540551240060986e-09, 4.639077530349239e-11, 9.446912374131907e-13, 1.239741169314320e-14, 9.821168268011667e-17, 4.281846844983401e-19, 8.943238836297170e-22, 7.129396795317248e-25, 1.410767108895213e-28, 2.344614766311755e-33};

static rl x24[24]={5.901985218150174e-02, 3.112391461984856e-01, 7.660969055459310e-01, 1.425597590803607e+00, 2.292562058632190e+00, 3.370774264209004e+00, 4.665083703467174e+00, 6.181535118736767e+00, 7.927539247172140e+00, 9.912098015077710e+00, 1.214610271172979e+01, 1.464273228959667e+01, 1.741799264650893e+01, 2.049146008261643e+01, 2.388732984816969e+01, 2.763593717433270e+01, 3.177604135237477e+01, 3.635840580165165e+01, 4.145172048487076e+01, 4.715310644515630e+01, 5.360857454469508e+01, 6.105853144721878e+01, 6.996224003510504e+01, 8.149827923394875e+01};
static rl w24[24]={1.428119733354791e-01, 2.587741075173471e-01, 2.588067072729614e-01, 1.833226889778158e-01, 9.816627262991873e-02, 4.073247815140497e-02, 1.322601940511963e-02, 3.369349058478236e-03, 6.721625640935950e-04, 1.044612146592732e-04, 1.254472197799206e-05, 1.151315812737281e-06, 7.960812959134761e-08, 4.072858987549941e-09, 1.507008226292724e-10, 3.917736515058580e-12, 6.894181052957633e-14, 7.819800382459149e-16, 5.350188813010072e-18, 2.010517464555542e-20, 3.605765864552929e-23, 2.451818845878373e-26, 4.088301593680631e-30, 5.575345788328810e-35};

static rl x25[25]={5.670477545269892e-02, 2.990108985869824e-01, 7.359095554350135e-01, 1.369183116035187e+00, 2.201326053721461e+00, 3.235675803558021e+00, 4.476496615073817e+00, 5.929083762700444e+00, 7.599899309956750e+00, 9.496749220932431e+00, 1.162901491177874e+01, 1.400795797654506e+01, 1.664712559728879e+01, 1.956289801146902e+01, 2.277524198683502e+01, 2.630877239096886e+01, 3.019429116331610e+01, 3.447109757192204e+01, 3.919060880393737e+01, 4.442234933616200e+01, 5.026457499383360e+01, 5.686496717394014e+01, 6.446667061595416e+01, 7.353423479210014e+01, 8.526015556249587e+01};
static rl w25[25]={1.375260142301543e-01, 2.516452737651651e-01, 2.561760028098040e-01, 1.862154903624783e-01, 1.031998481075360e-01, 4.471416112994486e-02, 1.530523288639853e-02, 4.152414632877233e-03, 8.920990732596833e-04, 1.511560191642437e-04, 2.006553180193398e-05, 2.067774396431962e-06, 1.634652022291098e-07, 9.766015062125436e-09, 4.327720794185161e-10, 1.389600963389560e-11, 3.138922792539945e-13, 4.802614822604183e-15, 4.735885364807537e-17, 2.814205379843084e-19, 9.164954395990587e-22, 1.418940009497273e-24, 8.273651944099053e-28, 1.168881711542775e-31, 1.315831500059196e-36};

static rl x26[26]={5.456448271791470e-02, 2.877079791061281e-01, 7.080160194786365e-01, 1.317081366043842e+00, 2.117120944628476e+00, 3.111093878877939e+00, 4.302770876040372e+00, 5.696818826169470e+00, 7.298909963070189e+00, 9.115862897294209e+00, 1.115582495493753e+01, 1.342850903245921e+01, 1.594550394814333e+01, 1.872068614546028e+01, 2.177077460854145e+01, 2.511609366231676e+01, 2.878164687345021e+01, 3.279867321840999e+01, 3.720698262961286e+01, 4.205861592443456e+01, 4.742389935628151e+01, 5.340218519725330e+01, 6.014277569007692e+01, 6.789147970545130e+01, 7.711799906933292e+01, 8.902840275041109e+01};
static rl w26[26]={1.326168841468234e-01, 2.448673680160829e-01, 2.533920221272899e-01, 1.886459836417798e-01, 1.079123416425698e-01, 4.865656544891800e-02, 1.748213797848933e-02, 5.024983836831110e-03, 1.155269160749565e-03, 2.117911583003684e-04, 3.079252604171158e-05, 3.523440091024355e-06, 3.141556325359097e-07, 2.155456814403754e-08, 1.120478117670427e-09, 4.329002405136713e-11, 1.213689159837112e-12, 2.396074223261228e-14, 3.205218559831475e-16, 2.762494150607840e-18, 1.433064364072763e-20, 4.064844694113018e-23, 5.459513322798759e-26, 2.743220987135306e-29, 3.300230478151928e-33, 3.083752269585364e-38};

static rl x27[27]={5.257989820635878e-02, 2.772291059877001e-01, 6.821639113780035e-01, 1.268814180750760e+00, 2.039159278250483e+00, 2.995835581875165e+00, 4.142194482563245e+00, 5.482371716049071e+00, 7.021376024605552e+00, 8.765202716850725e+00, 1.072097904599686e+01, 1.289715090190977e+01, 1.530372419706893e+01, 1.795258016401067e+01, 2.085789274347120e+01, 2.403669038598311e+01, 2.750962763367572e+01, 3.130207079202107e+01, 3.544567066613283e+01, 3.998072262112672e+01, 4.495986478533800e+01, 5.045419615272770e+01, 5.656413090151837e+01, 6.344054668468885e+01, 7.133184805516014e+01, 8.071276384709610e+01, 9.280261352555702e+01};
static rl w27[27]={1.280457269074737e-01, 2.384188493943396e-01, 2.504934065230675e-01, 1.906620675935009e-01, 1.123109493118593e-01, 5.253938547513114e-02, 1.974092327574478e-02, 5.984113555700753e-03, 1.464150360922784e-03, 2.884980039658172e-04, 4.557654440373148e-05, 5.735756119190609e-06, 5.702097017139201e-07, 4.430769157631981e-08, 2.656165316493600e-09, 1.208947524768718e-10, 4.096264427590561e-12, 1.008349767316482e-13, 1.749051751115433e-15, 2.055967083261440e-17, 1.556295736723773e-19, 7.081265177901326e-22, 1.757472010428563e-24, 2.057043036797539e-27, 8.947419115368966e-31, 9.209310916917627e-35, 7.179742378294809e-40};

static rl x28[28]={5.073462484987284e-02, 2.674872686407407e-01, 6.581366283547904e-01, 1.223971808384914e+00, 1.966767612473777e+00, 2.888883326030319e+00, 3.993311659250103e+00, 5.283736062843449e+00, 6.764603404243513e+00, 8.441216328271322e+00, 1.031985046299324e+01, 1.240790341446065e+01, 1.471408516413573e+01, 1.724866341560802e+01, 2.002378332995176e+01, 2.305389013503029e+01, 2.635629737440133e+01, 2.995196683359622e+01, 3.386660551658440e+01, 3.813225441019465e+01, 4.278967237077264e+01, 4.789207163362278e+01, 5.351129795966425e+01, 5.974879608464124e+01, 6.675697728390648e+01, 7.478677815233922e+01, 8.431783710722699e+01, 9.658242062752728e+01};
static rl w28[28]={1.237788439543970e-01, 2.322792769009197e-01, 2.475118960365057e-01, 1.923071131323180e-01, 1.164053617211324e-01, 5.634590536448111e-02, 2.206636432626232e-02, 7.025887635583455e-03, 1.820607892695724e-03, 3.833443038571280e-04, 6.535087080694968e-05, 8.971362053411803e-06, 9.847012256249947e-07, 8.564075852673943e-08, 5.836838763137529e-09, 3.075638877842322e-10, 1.232590952724392e-11, 3.682173674108047e-13, 7.998790575969710e-15, 1.224922500324068e-16, 1.271124295030564e-18, 8.488593367686220e-21, 3.402455379425650e-23, 7.420156588867466e-26, 7.600413205801689e-29, 2.873910317940254e-32, 2.541822903889406e-36, 1.661375878029070e-41};

static rl x29[29]={4.901448999471220e-02, 2.584072979187776e-01, 6.357472158006812e-01, 1.182201055624839e+00, 1.899366498002681e+00, 2.789363538356294e+00, 3.854876158195291e+00, 5.099200088698592e+00, 6.526302885192905e+00, 8.140899733223035e+00, 9.948548893634888e+00, 1.195577194404265e+01, 1.417020584332057e+01, 1.660079655231565e+01, 1.925804791313072e+01, 2.215434543487158e+01, 2.530438376305017e+01, 2.872574102526637e+01, 3.243966674232434e+01, 3.647218971363949e+01, 4.085572233672913e+01, 4.563146770127779e+01, 5.085319151599030e+01, 5.659346289666186e+01, 6.295472850114935e+01, 7.009089459642412e+01, 7.825537041941043e+01, 8.793259364298611e+01, 1.003674916027667e+02};
static rl w29[29]={1.197868672250461e-01, 2.264293889447646e-01, 2.444735653556417e-01, 1.936202644976772e-01, 1.202069996018019e-01, 6.006260605181780e-02, 2.444404002270963e-02, 8.145624870286855e-03, 2.225899096067055e-03, 4.982863213445047e-04, 9.111348999272606e-05, 1.354703463179014e-05, 1.627850005022406e-06, 1.568720699017842e-07, 1.200953310119946e-08, 7.220528067633611e-10, 3.362401315847788e-11, 1.192547905662338e-12, 3.156204723862947e-14, 6.078331294228438e-16, 8.254298347471802e-18, 7.594039930710428e-20, 4.492421504929952e-22, 1.592636814587942e-24, 3.064007640463017e-27, 2.757266953598470e-30, 9.099251850437396e-34, 6.943730710008932e-38, 3.822247128465026e-43};

static rl x30[30]={4.740718054080234e-02, 2.499239167531598e-01, 6.148334543927695e-01, 1.143195825666098e+00, 1.836454554622575e+00, 2.696521874557224e+00, 3.725814507779514e+00, 4.927293765849885e+00, 6.304515590965068e+00, 7.861693293370260e+00, 9.603775985479256e+00, 1.153654659795616e+01, 1.366674469306425e+01, 1.600222118898111e+01, 1.855213484014316e+01, 2.132720432178314e+01, 2.434003576453266e+01, 2.760555479678099e+01, 3.114158670111127e+01, 3.496965200824905e+01, 3.911608494906789e+01, 4.361365290848486e+01, 4.850398616380419e+01, 5.384138540650749e+01, 5.969912185923554e+01, 6.618061779443856e+01, 7.344123859555985e+01, 8.173681050672774e+01, 9.155646652253684e+01, 1.041575244310587e+02};
static rl w30[30]={1.160440860208239e-01, 2.208511247506685e-01, 2.413998275878685e-01, 1.946367684464437e-01, 1.237284159668692e-01, 6.367878036897413e-02, 2.686047527337868e-02, 9.338070881604063e-03, 2.680696891337074e-03, 6.351291219408808e-04, 1.239074599068909e-04, 1.982878843895139e-05, 2.589350929131293e-06, 2.740942840535626e-07, 2.332831165025712e-08, 1.580745574778349e-09, 8.427479123057731e-11, 3.485161234907757e-12, 1.099018059753409e-13, 2.588312664959322e-15, 4.437838059840301e-17, 5.365918308212121e-19, 4.393946892291795e-21, 2.311409794388686e-23, 7.274588498292221e-26, 1.239149701448193e-28, 9.832375083105892e-32, 2.842323553402687e-35, 1.878608031749574e-39, 8.745980440465938e-45};

static rl x31[31]={4.590194762110553e-02, 2.419801638247699e-01, 5.952538942223508e-01, 1.106689499533003e+00, 1.777595692874775e+00, 2.609703415256674e+00, 3.605196802340039e+00, 4.766747084471759e+00, 6.097554567181732e+00, 7.601400949233132e+00, 9.282714313470898e+00, 1.114664975561930e+01, 1.319918957624502e+01, 1.544726831554933e+01, 1.789892982664477e+01, 2.056352633671584e+01, 2.345197348201184e+01, 2.657708135211826e+01, 2.995399087234650e+01, 3.360075953290220e+01, 3.753916440733040e+01, 4.179583087018217e+01, 4.640386680641110e+01, 5.140531447679779e+01, 5.685499286871588e+01, 6.282685590878635e+01, 6.942527719108028e+01, 7.680704776386274e+01, 8.523035860754578e+01, 9.518893989152571e+01, 1.079522438275785e+02};
static rl w31[31]={1.125278955042064e-01, 2.155276081809917e-01, 2.383082516456806e-01, 1.953883092978481e-01, 1.269828328930492e-01, 6.718616892390888e-02, 2.930322499388210e-02, 1.059756991529595e-02, 3.185127258238988e-03, 7.954954830794413e-04, 1.648005212663575e-04, 2.822923786430920e-05, 3.980290255100360e-06, 4.593183984179813e-07, 4.307554518772849e-08, 3.255124993826993e-09, 1.962024667541173e-10, 9.319049908661848e-12, 3.437754181940755e-13, 9.679524713044675e-15, 2.036806611011679e-16, 3.121268728071578e-18, 3.372958170416228e-20, 2.467279638661563e-22, 1.158220190452523e-24, 3.247292259142455e-27, 4.914301730806028e-30, 3.450007110480820e-33, 8.766371011715471e-37, 5.036364392115892e-41, 1.990998458253363e-46};

static rl x32[32]={4.448936583326735e-02, 2.345261095196172e-01, 5.768846293018862e-01, 1.072448753817818e+00, 1.722408776444648e+00, 2.528336706425795e+00, 3.492213273021994e+00, 4.616456769749769e+00, 5.903958504174239e+00, 7.358126733186229e+00, 8.982940924212608e+00, 1.078301863253999e+01, 1.276369798674273e+01, 1.493113975552254e+01, 1.729245433671532e+01, 1.985586094033606e+01, 2.263088901319680e+01, 2.562863602245917e+01, 2.886210181632343e+01, 3.234662915396478e+01, 3.610049480575203e+01, 4.014571977153942e+01, 4.450920799575493e+01, 4.922439498730869e+01, 5.433372133339699e+01, 5.989250916213401e+01, 6.597537728793506e+01, 7.268762809066263e+01, 8.018744697791351e+01, 8.873534041789239e+01, 9.882954286828397e+01, 1.117513980979374e+02};
static rl w32[32]={1.092183419523226e-01, 2.104431079388359e-01, 2.352132296698536e-01, 1.959033359728866e-01, 1.299837862860609e-01, 7.057862386571742e-02, 3.176091250917516e-02, 1.191821483483833e-02, 3.738816294611733e-03, 9.808033066150542e-04, 2.148649188013465e-04, 3.920341967987604e-05, 5.934541612868495e-06, 7.416404578668029e-07, 7.604567879120734e-08, 6.350602226625737e-09, 4.281382971040667e-10, 2.305899491891820e-11, 9.799379288728057e-13, 3.237801657729005e-14, 8.171823443419926e-16, 1.542133833393881e-17, 2.119792290163628e-19, 2.054429673787901e-21, 1.346982586637255e-23, 5.661294130397396e-26, 1.418560545463026e-28, 1.913375494454364e-31, 1.192248760098234e-34, 2.671511219240155e-38, 1.338616942106261e-42, 4.510536193899647e-48};

static rl x33[33]={4.316113561733012e-02, 2.275178028033701e-01, 5.596166558515427e-01, 1.040268507751007e+00, 1.670559196075715e+00, 2.451920795897632e+00, 3.386155337588011e+00, 4.475459498399776e+00, 5.722454720272106e+00, 7.130224344400094e+00, 8.702359230621406e+00, 1.044301365020597e+01, 1.235697375935028e+01, 1.444974168158557e+01, 1.672763921863833e+01, 1.919793658721247e+01, 2.186901352492818e+01, 2.475056290615775e+01, 2.785385111141333e+01, 3.119205554557506e+01, 3.478070915353836e+01, 3.863829671777401e+01, 4.278707207825347e+01, 4.725420660299330e+01, 5.207345190151430e+01, 5.728763454109291e+01, 6.295256594690657e+01, 6.914351338010977e+01, 7.596668701424693e+01, 8.358163722327083e+01, 9.225113944413515e+01, 1.024778443368234e+02, 1.155547564489955e+02};
static rl w33[33]={1.060977455531555e-01, 2.055829836620170e-01, 2.321265234960045e-01, 1.962073727690839e-01, 1.327448567051734e-01, 7.385180388770883e-02, 3.422323341086037e-02, 1.329397518080761e-02, 4.340943095046144e-03, 1.192250990668848e-03, 2.751582255823986e-04, 5.324334099228301e-05, 8.609571326463825e-06, 1.158377961024547e-06, 1.289851148565202e-07, 1.180967869803168e-08, 8.822766409670507e-10, 5.329612134103736e-11, 2.575504037483360e-12, 9.831413322522208e-14, 2.920414955565337e-15, 6.630771567524034e-17, 1.126098637049974e-18, 1.393116571223850e-20, 1.214810098915323e-22, 7.161581811421279e-25, 2.703207124881349e-27, 6.071923612869965e-30, 7.321342111326654e-33, 4.061317061455856e-36, 8.049522845451767e-40, 3.529029903604569e-44, 1.017166562994320e-49};

static rl x34[34]={4.190992002006858e-02, 2.209164019523691e-01, 5.433536945565451e-01, 1.009967765484322e+00, 1.621751953473947e+00, 2.380014626242269e+00, 3.286400154001771e+00, 4.342910167051674e+00, 5.551929289036553e+00, 6.916256711692250e+00, 8.439144795688373e+00, 1.012434610894347e+01, 1.197617070083502e+01, 1.399955594701638e+01, 1.620015202849431e+01, 1.858442710740986e+01, 2.115979764987493e+01, 2.393479130658661e+01, 2.691925258122363e+01, 3.012460565236799e+01, 3.356419491650055e+01, 3.725373335110333e+01, 4.121190385773085e+01, 4.546118330644089e+01, 5.002900054129942e+01, 5.494941289395278e+01, 6.026562168468019e+01, 6.603391492353073e+01, 7.233019307779107e+01, 7.926155448779583e+01, 8.698888681494736e+01, 9.577719042218814e+01, 1.061334484823827e+02, 1.193621066777033e+02};
static rl w34[34]={1.031503857574030e-01, 2.009336243365965e-01, 2.290577133033010e-01, 1.963233093513149e-01, 1.352794683850790e-01, 7.700290130862526e-02, 3.668093396701032e-02, 1.471880421628392e-02, 4.990295392628834e-03, 1.430809946564837e-03, 3.467057685540947e-04, 7.087155216352575e-05, 1.218662032646028e-05, 1.756124027700453e-06, 2.110776099032411e-07, 2.104140883091538e-08, 1.727915984033995e-09, 1.159678783102055e-10, 6.301977516552469e-12, 2.742816924614809e-13, 9.438765904262911e-15, 2.529440327462427e-16, 5.183705189158635e-18, 7.947967859944002e-20, 8.876789005199131e-22, 6.985695546383411e-24, 3.713786539776405e-26, 1.262566747466280e-28, 2.549456478722615e-31, 2.755788341805534e-34, 1.364760308970231e-37, 2.399497948459281e-41, 9.232070476190538e-46, 2.283838051319149e-51};

static rl x35[35]={4.072920906171643e-02, 2.146874527351438e-01, 5.280103843193471e-01, 9.813861734590834e-01, 1.575725947577454e+00, 2.312228296513214e+00, 3.192397939493120e+00, 4.218064125029714e+00, 5.391402734673820e+00, 6.714963279915733e+00, 8.191701758837640e+00, 9.825020497092684e+00, 1.161881638890732e+01, 1.357753935730782e+01, 1.570626339576274e+01, 1.801077328783511e+01, 2.049767110714663e+01, 2.317450799779923e+01, 2.604994871037019e+01, 2.913397920954446e+01, 3.243817183767790e+01, 3.597602876989391e+01, 3.976343410480025e+01, 4.381926011860664e+01, 4.816619797401877e+01, 5.283192505801566e+01, 5.785079502213443e+01, 6.326637367538488e+01, 6.913541388364827e+01, 7.553443513472999e+01, 8.257140552358268e+01, 9.040852386440075e+01, 9.931297365944137e+01, 1.097959909449066e+02, 1.231732531753757e+02};
static rl w35[35]={1.003622374186903e-01, 1.964823833618064e-01, 2.260145669342552e-01, 1.962716678963576e-01, 1.376007420131165e-01, 8.003040046413989e-02, 3.912577142888370e-02, 1.618672652612652e-02, 5.685325571774636e-03, 1.697223299126086e-03, 4.304850681416093e-04, 9.263405306366882e-05, 1.687026486931568e-05, 2.591672147427802e-06, 3.344610957266900e-07, 3.607785438691756e-08, 3.233569032210266e-09, 2.391322981592309e-10, 1.447332399039025e-11, 7.101355277811058e-13, 2.793369509321816e-14, 8.694840214332735e-16, 2.108845936777249e-17, 3.912910811923565e-19, 5.432719406706777e-21, 5.493694724075208e-23, 3.912692565636496e-25, 1.880986928220020e-27, 5.775155104393900e-30, 1.051140746780529e-32, 1.021277013187373e-35, 4.527347312081734e-39, 7.080232733925808e-43, 2.397494633364857e-47, 5.106723813701156e-53};

static rl x36[36]={3.961320640860452e-02, 2.088002856543561e-01, 5.135107755351906e-01, 9.543811530791473e-01, 1.532249225393836e+00, 2.248215817854413e+00, 3.103661490498033e+00, 4.100262545699017e+00, 5.240010108754690e+00, 6.525233330467676e+00, 7.958627508998753e+00, 9.543288035278511e+00, 1.128275128666267e+01, 1.318104390024876e+01, 1.524274226656034e+01, 1.747304463141492e+01, 1.987785893307996e+01, 2.246391051439111e+01, 2.523887525730234e+01, 2.821154567688228e+01, 3.139204037421494e+01, 3.479207144800394e+01, 3.842529076480125e+01, 4.230774567207645e+01, 4.645849004290486e+01, 5.090042150332348e+01, 5.566145791468383e+01, 6.077624068452877e+01, 6.628869068025382e+01, 7.225601475445673e+01, 7.875533816278156e+01, 8.589548142949653e+01, 9.383992978484154e+01, 1.028580101471193e+02, 1.134651354897073e+02, 1.269880151966414e+02};
static rl w36[36]={9.772074843097506e-02, 1.922175131097200e-01, 2.230033446978749e-01, 1.960708469294170e-01, 1.397213896525591e-01, 8.293386566453573e-02, 4.155046205905646e-02, 1.769190808860058e-02, 6.424205254347524e-03, 1.992009755530736e-03, 5.274125390015615e-04, 1.190928293808261e-04, 2.288738819700135e-05, 3.732771855991790e-06, 5.147593556930683e-07, 5.975678522639870e-08, 5.808893359619386e-09, 4.699410135050238e-10, 3.141361420395762e-11, 1.720657768792050e-12, 7.648272252595379e-14, 2.727833331557431e-15, 7.703901826650590e-17, 1.696146717146739e-18, 2.857649674150513e-20, 3.602797339948777e-22, 3.307592946790374e-24, 2.137610078052896e-26, 9.316746849295045e-29, 2.589919452202921e-31, 4.259690662859414e-34, 3.729330417638632e-37, 1.483586406885461e-40, 2.069056568490571e-44, 6.182834840807875e-49, 1.137392833407645e-54};

static rl x37[37]={3.855673418601471e-02, 2.032275099409546e-01, 4.997870664422638e-01, 9.288254993469235e-01, 1.491115012042621e+00, 2.187669071935328e+00, 3.019757467642508e+00, 3.988920302941583e+00, 5.096984544604990e+00, 6.346084060649569e+00, 7.738683813536452e+00, 9.277608941806196e+00, 1.096607963131078e+01, 1.280775291254541e+01, 1.480677283245775e+01, 1.696783086160973e+01, 1.929623894818940e+01, 2.179801837446761e+01, 2.448000859074969e+01, 2.735000161963780e+01, 3.041690962745649e+01, 3.369097614342631e+01, 3.718404563922262e+01, 4.090991252978045e+01, 4.488478041515481e+01, 4.912787780671970e+01, 5.366230173869413e+01, 5.851620322014998e+01, 6.372450357544494e+01, 6.933147002696597e+01, 7.539475330169681e+01, 8.199207200019234e+01, 8.923308272823333e+01, 9.728253366634951e+01, 1.064118573194836e+02, 1.171405711133349e+02, 1.308062253516501e+02};
static rl w37[37]={9.521465142709475e-02, 1.881281008672653e-01, 2.200290515515567e-01, 1.957373423785228e-01, 1.416536425236323e-01, 8.571375643619764e-02, 4.394862142284201e-02, 1.922871087836737e-02, 7.204877170310809e-03, 2.315471642949444e-03, 6.383325586365018e-04, 1.508182894937913e-04, 3.048561891686146e-05, 5.258684995135365e-06, 7.716106400103869e-07, 9.592841294453468e-08, 1.005754196237682e-08, 8.844088579626155e-10, 6.481397346488567e-11, 3.929547811692776e-12, 1.954245790121780e-13, 7.894002371087710e-15, 2.560463596349619e-16, 6.579947209691782e-18, 1.318719306639603e-19, 2.022821053234999e-21, 2.321861994421940e-23, 1.940161521429577e-25, 1.140611623635595e-27, 4.518102148685397e-30, 1.139888100662706e-32, 1.698157942248677e-35, 1.342849603575171e-38, 4.805349653642180e-42, 5.990990619381185e-46, 1.583922443954436e-50, 2.523773190544875e-56};

static rl x38[38]={3.755515263135437e-02, 1.979445866473918e-01, 4.867785379490556e-01, 9.046053592858760e-01, 1.452138374380256e+00, 2.130312744494631e+00, 2.940299106412451e+00, 3.883515862381367e+00, 4.961643594510520e+00, 6.176642448847499e+00, 7.530773048708940e+00, 9.026621537140498e+00, 1.066712947636635e+01, 1.245562951416466e+01, 1.439588796713035e+01, 1.649215578524130e+01, 1.874922977721838e+01, 2.117252653020712e+01, 2.376817220623610e+01, 2.654311242428542e+01, 2.950524786450047e+01, 3.266360324951338e+01, 3.602854026261596e+01, 3.961202922587127e+01, 4.342800075688229e+01, 4.749280845538381e+01, 5.182584921134675e+01, 5.645041307625637e+01, 6.139487752088784e+01, 6.669443657755892e+01, 7.239369575777199e+01, 7.855074022042079e+01, 8.524387022996996e+01, 9.258356300923502e+01, 1.007358073835304e+02, 1.099741053589464e+02, 1.208220095175347e+02, 1.346277282875351e+02};
static rl w38[38]={9.283380362872036e-02, 1.842040072834431e-01, 2.170956462404759e-01, 1.952859469798631e-01, 1.434092042396506e-01, 8.837126752057826e-02, 4.631470042464929e-02, 2.079173405822041e-02, 8.025103445714577e-03, 2.667706176052068e-03, 7.640088976752587e-04, 1.883818369473179e-04, 3.993130886529885e-05, 7.260522350193317e-06, 1.129175089417986e-06, 1.496825625111701e-07, 1.684167542133368e-08, 1.600600390696633e-09, 1.277633557456908e-10, 8.509826002853740e-12, 4.694201426342415e-13, 2.126021701262418e-14, 7.826939565222512e-16, 2.315208340320794e-17, 5.428404841374320e-19, 9.929073707616101e-21, 1.390171625266857e-22, 1.456330948053784e-24, 1.110270426662952e-26, 5.951511121348533e-29, 2.147456054081499e-31, 4.928296799396074e-34, 6.665155907588076e-37, 4.771205396858213e-40, 1.539295212785391e-43, 1.719550191139237e-47, 4.032092947462518e-52, 5.580031520358049e-58};

static rl x39[39]={3.660429195660375e-02, 1.929294667748622e-01, 4.744306514626067e-01, 8.816185217504263e-01, 1.415153402964431e+00, 2.075900054441475e+00, 2.864940087124226e+00, 3.783582804450770e+00, 4.833377810488478e+00, 6.016130155727397e+00, 7.333918496421592e+00, 8.789116761975695e+00, 1.038441992269649e+01, 1.212287459214670e+01, 1.400791538476811e+01, 1.604340819160565e+01, 1.823370185334919e+01, 2.058369012744060e+01, 2.309888640346517e+01, 2.578551437483859e+01, 2.865061890810396e+01, 3.170220278874495e+01, 3.494939705509712e+01, 3.840267555752546e+01, 4.207412867423611e+01, 4.597781755745090e+01, 5.013024018849510e+01, 5.455095617610832e+01, 5.926344277150436e+01, 6.429629777764541e+01, 6.968498127741620e+01, 7.547442956780553e+01, 8.172315339323799e+01, 8.851002366778249e+01, 9.594632356854197e+01, 1.041992613232035e+02, 1.135443739641761e+02, 1.245091828960139e+02, 1.384523795310612e+02};
static rl w39[39]={9.056905005139866e-02, 1.804358080434677e-01, 2.142062151332155e-01, 1.947299294910617e-01, 1.449992236200189e-01, 9.090819105555371e-02, 4.864391970979148e-02, 2.237584372362728e-02, 8.882509779705101e-03, 3.048619104381225e-03, 9.051184193678579e-04, 2.323487733734230e-04, 5.150713378858222e-05, 9.841350607160249e-06, 1.616603654561842e-06, 2.275914932550321e-07, 2.735843835983178e-08, 2.795839409850404e-09, 2.416700259571923e-10, 1.756695502203115e-11, 1.066680310239033e-12, 5.369209358741653e-14, 2.220765375245781e-15, 7.471439791559861e-17, 2.020732421263994e-18, 4.333726724871730e-20, 7.252087338658738e-22, 9.289780459158246e-24, 8.902481021644814e-26, 6.206183530706068e-28, 3.040039486340603e-30, 1.001373894925193e-32, 2.094905241211365e-35, 2.577494120748580e-38, 1.673809533334302e-41, 4.878962035553728e-45, 4.894358450929249e-49, 1.020239070240378e-53, 1.229529448418639e-59};

 static rl* xlookup[]={x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39};//, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50};

static rl* wlookup[]={w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17, w18, w19, w20, w21, w22, w23, w24, w25, w26, w27, w28, w29, w30, w31, w32, w33, w34, w35, w36, w37, w38, w39};//, w40, w41, w42, w43, w44, w45, w46, w47, w48, w49, w50};

int nmax=39;
if((n>0) and (n<=nmax)){
  for(int i=0;i<n;i++){
    glx[i]=xlookup[n-1][i];
    glw[i]=wlookup[n-1][i];
 }
  return;
 }
 cout<< "n out of range in laguerre_table!! n="<<n<<"> nmax="<<nmax<<"\n";
}