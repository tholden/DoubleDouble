% The below code is derived from the QD C++ library.

% QD is Copyright ( c ) 2003-2009, The Regents of the University of California, through Lawrence Berkeley National Laboratory ( subject to receipt of any required approvals from U.S. Dept. of Energy ) All rights reserved.

% QD is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% ( 3 ) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 3. You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ( "Enhancements" ) to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.

% The implementation of the LUP decomposition and the backslack operator here is derived from Cleve Moler's code from "Numerical Computing with MATLAB".

% "Numerical Computing with MATLAB" is Copyright ( c ) 2004, Cleve Moler and Copyright ( c ) 2016, The MathWorks, Inc.

% "Numerical Computing with MATLAB" is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% ( 3 ) In all cases, the software is, and all modifications and derivatives of the software shall be, licensed to you solely for use in conjunction with MathWorks products and service offerings.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% The implementation of the LDL decomposition here is derived from Brian Borchers's ldlt library.

% Brian Borchers's ldlt library is Copyright ( c ) 2009, Brian Borchers.

% Brian Borchers's ldlt library is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
% ( 3 ) Neither the name of the New Mexico Inst of Mining & Tech nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% WARNING: This code was primarily written by AI, based on OctDouble.m, using DoubleDouble components. Having not gone through it line by line, I cannot 100% guarantee its correctness.

classdef QuadDoubleSlow < BaseExtDouble

    properties ( Constant, GetAccess = public )

        empty  = QuadDoubleSlow.MakeStatic( DoubleDouble.empty, DoubleDouble.empty );
        zero   = QuadDoubleSlow.MakeStatic( DoubleDouble.zero, DoubleDouble.zero );
        one    = QuadDoubleSlow.MakeStatic( DoubleDouble.one, DoubleDouble.zero );
        tiny   = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.21543267145725000000e-63, 0.00000000000000000000e+00 ), DoubleDouble.zero );
        pi     = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.14159265358979311600e+00, 1.22464679914735320717e-16 ), DoubleDouble.MakeStatic( -2.99476980971833966589e-33, 1.11245422086336528166e-49 ) );
        piT2   = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 6.28318530717958623200e+00, 2.44929359829470641444e-16 ), DoubleDouble.MakeStatic( -5.98953961943667933152e-33, 2.22490844172673056346e-49 ) );
        piD2   = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.57079632679489655800e+00, 6.12323399573676603611e-17 ), DoubleDouble.MakeStatic( -1.49738490485916983288e-33, 5.56227110431682640865e-50 ) );
        piD16  = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.96349540849362069750e-01, 7.65404249467095754484e-18 ), DoubleDouble.MakeStatic( -1.87173113107396229118e-34, 6.95283888039603300966e-51 ) );
        log_2  = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 6.93147180559945286227e-01, 2.31904681384629955842e-17 ), DoubleDouble.MakeStatic( 5.70770843841621206578e-34, -3.58243221060181142336e-50 ) );
        log_10 = QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.30258509299404590109e+00, -2.17075622338224935076e-16 ), DoubleDouble.MakeStatic( -9.98426245446577657012e-33, -4.02335745445020637879e-49 ) );

        NInverseFactorial = 48;
        ExpRescale        = 12;
        LogSteps          = 3;

        InverseFactorial = [
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.66666666666666657415e-01, 9.25185853854297065662e-18 ), DoubleDouble.MakeStatic( 5.13581318503262865639e-34, 2.85094902409834186429e-50 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 4.16666666666666643537e-02, 2.31296463463574266415e-18 ), DoubleDouble.MakeStatic( 1.28395329625815716410e-34, 7.12737256024585466073e-51 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.33333333333333321769e-03, 1.15648231731787138023e-19 ), DoubleDouble.MakeStatic( 1.60494162032269652194e-36, 2.22730392507682967421e-53 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.38888888888888894189e-03, -5.30054395437357705906e-20 ), DoubleDouble.MakeStatic( -1.73868675534958775956e-36, -1.63335621172300839684e-52 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.98412698412698412526e-04, 1.72095582934207052868e-22 ), DoubleDouble.MakeStatic( 1.49269123913941270724e-40, 1.29470326746002470656e-58 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.48015873015873015658e-05, 2.15119478667758816085e-23 ), DoubleDouble.MakeStatic( 1.86586404892426588405e-41, 1.61837908432503088320e-59 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.75573192239858925110e-06, -1.85839327404647208104e-22 ), DoubleDouble.MakeStatic( 8.49175460488199287009e-39, -5.72661640789429621316e-55 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.75573192239858882758e-07, 2.37677146222502973185e-23 ), DoubleDouble.MakeStatic( -3.26318890334088294370e-40, 1.61435111860404415106e-56 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.50521083854417202239e-08, -1.44881407093591196603e-24 ), DoubleDouble.MakeStatic( 2.04267351467144545891e-41, -8.49632672007163174711e-58 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.08767569878681001866e-09, -1.20734505911325997169e-25 ), DoubleDouble.MakeStatic( 1.70222792889287100335e-42, 1.41609532150396699816e-58 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.60590438368216133409e-10, 1.25852945887520980521e-26 ), DoubleDouble.MakeStatic( -5.31334602762985030694e-43, 3.54021472597605527826e-59 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.14707455977297245073e-11, 2.06555127528307454245e-28 ), DoubleDouble.MakeStatic( 6.88907923246664603290e-45, 5.72920002655109095474e-61 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 7.64716373181981640551e-13, 7.03872877733453001061e-30 ), DoubleDouble.MakeStatic( -7.82753927716258344520e-48, 1.92138649443790241643e-64 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 4.77947733238738525345e-14, 4.39920548583408125663e-31 ), DoubleDouble.MakeStatic( -4.89221204822661465325e-49, 1.20086655902368901027e-65 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.81145725434552059811e-15, 1.65088427308614325994e-31 ), DoubleDouble.MakeStatic( -2.87777179307447917987e-50, 4.27110689256293549037e-67 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.56192069685862252711e-16, 1.19106796602737540024e-32 ), DoubleDouble.MakeStatic( -4.57750605962998323416e-49, 2.87494142340899603160e-67 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.22063524662432949554e-18, 2.21418941196042653637e-34 ), DoubleDouble.MakeStatic( -1.50891402377419897072e-50, 1.40072951514781547649e-67 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 4.11031762331216484407e-19, 1.44129733786595271498e-36 ), DoubleDouble.MakeStatic( -5.28562754878981208303e-53, -4.14764725635765684990e-70 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.95729410633912625952e-20, -1.36435038300879084872e-36 ), DoubleDouble.MakeStatic( 1.33923482511250642308e-53, -6.82108942414933121893e-70 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.89679139245057407789e-22, -7.91140261487237621703e-38 ), DoubleDouble.MakeStatic( -3.18779767905709332675e-54, 1.27057810175205661633e-70 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.86817017063068412615e-23, -8.84317765548234384789e-40 ), DoubleDouble.MakeStatic( 3.87181571061732467175e-56, -1.95652575315225570181e-72 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.61173757109611838590e-24, -3.68465735645097660329e-41 ), DoubleDouble.MakeStatic( 1.61325654609055194656e-57, -8.15219063813439928119e-74 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 6.44695028438447358950e-26, -1.93304042337034648245e-42 ), DoubleDouble.MakeStatic( -1.52130238070391441718e-58, 6.64377273721295752858e-75 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.47959626322479758961e-27, -1.29537309647652287758e-43 ), DoubleDouble.MakeStatic( 6.40339015984996240505e-60, -8.46024562770674585087e-77 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 9.18368986379554600539e-29, 1.43031503967873220242e-45 ), DoubleDouble.MakeStatic( -8.55122677465050479948e-62, 8.38146710023453831785e-78 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.27988923706983775898e-30, 1.51175427440298786897e-46 ), DoubleDouble.MakeStatic( 8.05851771951971592852e-63, -9.09648053071092885400e-81 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.13099628864477158818e-31, 1.04980154129595060239e-47 ), DoubleDouble.MakeStatic( -4.34615092939779517622e-64, -4.96677980014005581500e-81 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.76998762881590538515e-33, 2.58703478327503238348e-49 ), DoubleDouble.MakeStatic( 3.23789002742563998619e-66, 2.56128591057885727338e-82 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.21612504155351789377e-34, 5.58629056788880576964e-51 ), DoubleDouble.MakeStatic( 6.61594857808279192552e-68, -3.16204422895208590789e-84 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.80039075485474341802e-36, 1.74571580246525180301e-52 ), DoubleDouble.MakeStatic( 2.06748393065087247673e-69, -9.88138821547526846215e-86 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.15163356207719508905e-37, -6.09957445788453977976e-54 ), DoubleDouble.MakeStatic( -5.34474961965941048489e-70, 2.62531262385000815484e-86 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.38715753552116179633e-39, 5.09056148151084994757e-56 ), DoubleDouble.MakeStatic( 3.98956734903634402553e-72, -1.14951294479092623352e-88 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 9.67759295863189067186e-41, 3.20229554864556196033e-57 ), DoubleDouble.MakeStatic( 6.54750720501810103533e-74, -5.91334284153607618996e-91 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.68822026628663633314e-42, 5.35506116594333401341e-59 ), DoubleDouble.MakeStatic( -1.12906019874498675666e-75, -7.09714352853527273815e-92 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 7.26546017915307135902e-44, -4.36409714935444569150e-61 ), DoubleDouble.MakeStatic( 2.55032501210183748394e-77, 3.62259693228430959850e-94 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.91196320504028195296e-45, -2.78608221768831260969e-62 ), DoubleDouble.MakeStatic( 2.03474372241013280071e-78, -9.13939362246162656853e-95 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 4.90246975651354351900e-47, -1.21301910051792795265e-63 ), DoubleDouble.MakeStatic( -4.47071800113765855194e-80, 5.37597340717858998037e-97 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.22561743912838584936e-48, 6.03392734831560538813e-68 ), DoubleDouble.MakeStatic( 4.07624961245823729659e-84, -1.59820353307989972028e-102 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 2.98931082714240461486e-50, -1.04072477030331555467e-66 ), DoubleDouble.MakeStatic( -1.37613197137759057242e-83, -5.01838323088147032117e-100 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 7.11740673129143899398e-52, 3.17420753842055730074e-68 ), DoubleDouble.MakeStatic( 1.24112898646225877483e-84, -1.09918879326293740961e-100 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.65521086774219514723e-53, 4.14710519049482419382e-70 ), DoubleDouble.MakeStatic( 4.32187417758181305900e-88, 1.61954498080637493977e-104 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.76184288123226156652e-55, 2.25971359112361835283e-71 ), DoubleDouble.MakeStatic( -1.45255187387709499120e-87, -1.23159968713973768582e-104 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.35965084718280448735e-57, -5.04027988508830641661e-73 ), DoubleDouble.MakeStatic( -1.97116512594399310946e-89, -1.18661232656178290706e-106 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.81731540156147899030e-58, 1.36506933987936602796e-74 ), DoubleDouble.MakeStatic( 2.54490150401424344048e-91, -2.39064169520447011495e-107 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.86662851396059404326e-60, -1.56435500578638897666e-76 ), DoubleDouble.MakeStatic( -1.28638554473548469031e-92, -4.72409213107264703206e-109 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.05547607075123643865e-62, 8.25581847807094913248e-78 ), DoubleDouble.MakeStatic( -2.67996988486559310482e-94, -9.84185860640134798346e-111 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.64397470831657907432e-63, -4.08088098184429380620e-80 ), DoubleDouble.MakeStatic( -3.32913139064192010208e-96, 1.47073743046046679366e-112 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.28794941663315801371e-65, 5.33225140364648136571e-82 ), DoubleDouble.MakeStatic( 8.32419386223677346691e-99, 6.67478974686111897684e-115 ) );
            ];

        SinTable = [
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 0.00000000000000000000e+00, 0.00000000000000000000e+00 ), DoubleDouble.MakeStatic( 0.00000000000000000000e+00, 0.00000000000000000000e+00 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.95090322016128275839e-01, -7.99107906846173126344e-18 ), DoubleDouble.MakeStatic( 6.18462700242207127059e-34, -3.58402709180329369584e-50 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 3.82683432365089781779e-01, -1.00507726964615876117e-17 ), DoubleDouble.MakeStatic( -2.06053163028066946701e-34, -1.27177246980852050275e-50 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 5.55570233019602177649e-01, 4.70941094056167682138e-17 ), DoubleDouble.MakeStatic( -2.06405203836829206178e-33, 1.22901631885671376629e-49 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 7.07106781186547572737e-01, -4.83364665672645672553e-17 ), DoubleDouble.MakeStatic( 2.06933765434970678136e-33, 2.46777349573417545616e-50 ) );
            ];

        CosTable = [
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 1.00000000000000000000e+00, 0.00000000000000000000e+00 ), DoubleDouble.MakeStatic( 0.00000000000000000000e+00, 0.00000000000000000000e+00 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 9.80785280403230430579e-01, 1.85469399978250057259e-17 ), DoubleDouble.MakeStatic( -1.06965644455307566037e-33, 6.66681744752649605784e-50 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 9.23879532511286738483e-01, 1.76450470843366770600e-17 ), DoubleDouble.MakeStatic( -5.04425373215868178317e-34, -4.04786777168238900471e-50 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 8.31469612302545235671e-01, 1.40738569847280238931e-18 ), DoubleDouble.MakeStatic( 4.69513153839808352459e-35, -2.02338815193825684523e-52 ) );
            QuadDoubleSlow.MakeStatic( DoubleDouble.MakeStatic( 7.07106781186547572737e-01, -4.83364665672645672553e-17 ), DoubleDouble.MakeStatic( 2.06933765434970678136e-33, 2.46777349573417545616e-50 ) );
            ];

    end

    methods

        function v = QuadDoubleSlow( in, varargin )
            if nargin == 0
                v = QuadDoubleSlow.empty;
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = QuadDoubleSlow.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDoubleSlow' ) || isa( in, 'QuadDouble' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            elseif isa( in, 'DoubleDouble' )
                v.v1 = in;
                v.v2 = DoubleDouble.zeros( size( in ) );
            elseif isa( in, 'BaseExtDouble' )
                C = cell( 1, 4 );
                [ C{ : } ] = ToSumOfDoubles( in );
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );
            else
                v.v1 = DoubleDouble( in );
                v.v2 = DoubleDouble.zeros( size( in ) );
            end
        end

        function v = Promote( ~, v )
            v = QuadDoubleSlow( v );
        end

        function n = PromotionOrder( ~ )
            n = 2;
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDoubleSlow.MakeStatic( a1, a2 );
        end

    end

    methods ( Static )

        function v = ones( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( ones( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = zeros( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( zeros( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = eye( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( eye( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = NaN( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble.NaN( varargin{:} ), DoubleDouble.NaN( varargin{:} ) );
        end

        function v = Inf( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble.Inf( varargin{:} ), DoubleDouble.Inf( varargin{:} ) );
        end

        function v = randn( varargin )
            v = ToRandn( QuadDoubleSlow.zeros( varargin{:} ) );
        end

        function v = rand( varargin )
            v = ToRand( QuadDoubleSlow.zeros( varargin{:} ) );
        end

        function v = randi( imax, varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( randi( imax, varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

    end

    methods ( Static, Access = { ?QuadDoubleSlow, ?QuadDouble, ?OctDouble } )

        function v = MakeStatic( a1, a2 )
            v = QuadDoubleSlow;
            if ~isa( a1, 'DoubleDouble' )
                a1 = DoubleDouble( a1 );
            end
            if ~isa( a2, 'DoubleDouble' )
                a2 = DoubleDouble( a2 );
            end
            v.v1 = a1;
            v.v2 = a2;
        end

    end

end
