﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FITS_READER
{
    class ThAr
    {
        public static double[] Wls = new double[]{
             3184.9486,3188.2339,3204.3210,3205.0912,3210.3085,3211.1943,3212.0220,
             3213.5741,3220.3508,3227.8028,3229.0096,3234.9962,3235.8400,3238.9343,
             3241.1079,3242.3790,3243.6887,3244.4488,3245.7603,3249.8003,3251.9159,
             3253.8659,3256.2738,3257.3667,3258.1220,3262.6684,3272.0268,3273.9157,
             3275.0676,3276.1712,3280.3713,3281.7016,3282.6166,3285.7525,3286.5829,
             3287.7893,3291.7394,3292.5209,3293.6988,3298.0498,3299.6687,3304.2383,
             3305.3036,3307.2283,3309.3654,3313.6784,3314.8268,3321.4508,3322.0933,
             3327.1931,3329.7284,3330.4770,3333.1290,3334.6041,3336.1618,3337.8703,
             3340.7254,3346.5560,3348.7684,3354.1796,3358.6020,3360.9982,3364.6855,
             3365.3383,3366.5171,3371.7967,3372.8230,3373.4925,3374.9749,3376.4359,
             3378.5734,3380.8595,3383.1068,3385.5316,3386.5006,3388.5309,3392.0349,
             3396.7278,3397.5161,3398.5448,3401.7110,3402.6952,3405.5584,3408.7499,
             3410.0756,3413.0130,3415.8846,3417.4978,3421.2100,3422.6561,3423.9897,
             3427.0923,3428.9992,3431.8104,3433.9988,3434.7271,3435.9771,3436.7272,
             3437.3071,3439.3987,3442.5790,3446.5474,3451.7023,3452.6820,3454.0952,
             3457.0691,3462.8505,3464.1272,3465.0626,3465.7650,3468.2198,3469.3454,
             3469.9208,3471.2186,3471.9593,3476.7474,3478.2324,3479.1725,3480.5055,
             3482.7613,3486.5512,3489.5076,3491.5360,3493.5185,3495.6998,3496.8107,
             3498.0098,3498.6210,3501.8666,3505.4946,3509.7785,3511.1574,3514.3877,
             3518.4040,3519.9936,3521.9139,3523.5061,3526.6342,3528.4116,3528.9544,
             3529.3859,3530.5148,3531.4505,3533.1826,3535.3196,3536.0108,3539.5872,
             3542.4979,3544.0179,3547.3376,3548.5144,3549.5959,3551.4019,3554.3058,
             3555.0135,3556.9041,3559.5081,3561.0304,3561.7809,3565.0298,3569.8204,
             3572.3923,3576.6156,3581.6084,3582.3546,3583.1022,3584.1756,3588.4407,
             3591.4524,3592.7794,3598.1199,3601.0344,3604.6819,3606.5218,3608.3770,
             3609.4452,3612.4275,3612.8666,3615.1327,3615.8503,3618.3633,3619.2119,
             3622.7954,3632.8303,3634.5822,3635.9433,3642.2490,3643.5123,3649.2496,
             3649.7349,3652.1683,3654.4618,3656.6938,3658.8087,3659.6294,3660.4370,
             3661.6214,3663.2025,3666.9811,3668.1398,3671.5398,3672.3003,3672.5218,
             3673.7935,3674.8910,3675.1372,3675.5675,3675.7893,3675.9592,3676.1424,
             3678.2701,3678.4804,3679.1343,3679.5444,3679.7105,3680.0609,3680.4477,
             3680.6053,3681.8836,3682.4863,3684.9329,3687.6700,3687.9841,3688.7604,
             3690.6238,3691.4117,3691.6141,3692.5664,3694.1785,3695.2889,3695.9737,
             3697.0308,3697.7436,3698.1061,3699.8808,3700.9872,3703.2299,3703.7743,
             3704.8616,3706.7672,3711.3041,3711.6229,3715.5606,3715.8615,3716.5836,
             3717.1713,3717.8293,3718.2065,3719.4347,3719.8362,3720.3073,3721.2152,
             3721.8254,3723.6561,3725.3932,3726.7246,3727.6120,3727.9027,3729.3087,
             3730.3683,3730.7485,3732.9854,3733.6724,3737.5125,3737.8890,3740.8551,
             3741.1830,3742.9234,3743.5083,3744.7367,3745.6591,3745.9706,3747.5390,
             3748.2838,3749.0843,3750.4939,3751.0219,3752.5700,3753.2421,3753.5177,
             3754.0308,3754.5930,3755.2121,3756.2941,3757.6941,3758.4671,3758.7063,
             3763.5053,3765.2700,3766.1186,3766.4473,3767.9007,3770.0560,3770.5200,
             3771.3708,3772.6498,3773.7574,3776.2711,3777.4167,3781.3189,3783.0127,
             3783.2964,3784.5749,3784.7942,3785.2804,3785.6002,3786.3824,3786.8830,
             3789.1680,3790.3560,3790.7950,3792.3740,3792.7300,3794.1510,3795.3800,
             3798.1030,3799.3820,3800.1980,3801.4430,3803.0750,3803.9840,3805.8200,
             3807.8740,3808.1290,3809.4561,3809.8350,3813.0680,3816.1666,3818.6860,
             3820.7926,3821.4308,3822.1485,3822.8619,3823.0676,3823.4577,3825.1331,
             3825.6729,3826.3688,3826.8072,3828.3846,3830.0606,3830.7736,3831.6398,
             3833.0860,3834.6787,3835.7111,3836.5851,3837.8752,3839.6953,3840.8004,
             3841.5187,3841.9601,3842.8968,3844.7311,3845.4055,3846.8876,3850.5813,
             3852.1353,3854.5108,3856.3544,3859.8396,3861.3524,3862.0522,3862.4218,
             3863.4059,3866.9092,3867.5786,3868.2547,3868.5284,3869.6633,3871.6174,
             3872.1371,3872.7222,3873.1479,3873.4740,3873.8224,3874.2438,3874.8619,
             3875.3731,3875.6462,3877.4630,3878.6620,3879.2688,3879.4454,3879.6441,
             3880.1948,3880.3332,3881.4980,3882.1430,3884.5245,3884.8225,3885.2254,
             3885.7681,3886.9159,3889.9059,3891.4017,3891.7260,3891.9792,3893.6519,
             3894.6600,3895.4192,3898.4374,3898.7967,3899.0310,3900.5767,3900.8785,
             3901.6621,3903.1024,3903.4812,3904.0828,3905.1865,3907.1603,3907.5432,
             3908.7491,3908.9797,3909.1390,3910.7734,3911.5760,3911.9091,3913.6450,
             3914.1623,3915.2953,3915.8489,3916.4176,3917.2693,3918.0701,3918.5115,
             3919.0234,3921.7799,3922.2191,3923.7995,3924.4034,3925.0934,3925.7188,
             3926.0456,3927.8062,3928.6233,3929.2904,3929.6693,3930.3303,3931.2566,
             3932.2264,3932.9113,3933.2375,3933.6611,3934.2742,3936.3532,3936.6718,
             3937.0406,3937.9237,3938.6144,3942.0732,3944.2534,3945.5072,3945.8201,
             3946.1455,3947.3308,3947.5046,3948.0305,3948.9789,3950.3951,3950.8050,
             3951.5149,3952.7608,3955.1699,3955.8903,3956.6908,3959.3000,3960.2692,
             3962.4196,3963.4688,3964.0302,3964.7257,3966.9647,3967.3921,3968.4673,
             3969.0026,3969.6641,3972.1545,3972.6395,3973.1961,3974.4766,3974.7590,
             3975.4674,3976.4149,3979.0433,3979.3559,3980.0896,3980.7550,3981.1065,
             3981.8281,3982.1020,3982.8956,3984.8795,3986.0739,3987.2061,3988.5997,
             3988.8456,3990.4922,3991.7309,3992.0535,3992.2726,3993.2989,3994.5494,
             3994.7918,3996.0617,3996.6687,3997.8654,3998.4038,3998.7334,3998.9526,
             4000.2811,4001.0581,4001.8935,4003.1057,4003.3086,4003.5740,4005.0928,
             4005.3628,4005.9615,4006.3806,4007.0188,4008.2102,4009.0112,4011.5916,
             4011.7398,4012.4952,4013.8566,4014.5143,4014.7164,4018.0990,4019.1289,
             4020.3541,4021.1498,4021.7506,4022.0674,4023.3378,4024.8025,4025.6556,
             4027.0091,4029.8256,4030.2925,4030.8424,4034.2461,4034.9218,4035.4600,
             4036.0479,4036.5652,4037.5614,4038.2287,4038.8043,4039.8648,4040.9155,
             4041.2036,4042.8937,4043.1302,4043.3948,4044.4179,4044.9269,4045.2268,
             4045.9654,4047.4815,4048.2876,4049.5361,4049.9447,4050.8872,4051.4991,
             4052.9208,4053.5277,4054.3017,4054.5258,4057.9412,4059.2529,4063.4071,
             4064.3315,4065.6187,4067.4507,4069.2014,4069.2833,4069.4612,4070.2383,
             4071.7513,4071.9967,4072.3849,4072.6284,4073.8563,4074.6420,4075.5030,
             4075.9070,4076.6284,4076.9432,4078.3015,4078.8751,4079.5738,4080.3580,
             4080.7058,4081.3678,4081.5914,4082.0817,4082.3872,4083.4688,4084.6192,
             4085.0421,4085.2564,4085.4341,4086.5205,4087.2848,4088.7264,4089.1379,
             4091.3476,4093.3921,4093.6713,4094.7470,4096.0759,4097.7478,4098.7316,
             4098.9327,4100.3414,4102.6177,4103.9121,4104.3821,4104.8385,4105.9150,
             4107.0506,4107.8614,4108.4198,4109.3234,4110.8260,4112.7545,4115.7589,
             4116.7137,4118.4893,4124.0403,4127.4120,4128.0499,4128.6400,4129.5116,
             4131.0021,4131.7235,4132.7533,4134.0681,4134.3235,4135.4800,4138.0406,
             4140.2350,4140.3855,4141.6321,4142.4741,4142.7010,4143.6491,4147.5396,
             4148.1816,4148.7269,4151.1453,4154.7205,4156.5166,4158.5905,4159.6599,
             4161.7389,4162.5090,4163.6469,4163.9479,4164.1795,4164.2550,4165.7661,
             4168.0464,4168.6341,4170.5333,4171.3410,4172.6257,4173.4475,4174.4396,
             4178.0597,4178.3658,4178.8481,4179.9601,4181.8836,4182.9342,4184.1376,
             4185.1447,4187.1410,4190.7129,4191.0294,4192.3620,4193.0164,4194.9318,
             4198.3170,4200.6745,4207.6122,4208.4108,4208.8907,4209.1964,4210.9232,
             4213.0673,4213.3779,4214.8285,4215.6272,4216.0698,4216.3722,4217.4308,
             4218.6649,4220.0651,4220.7323,4221.1530,4221.6921,4222.6373,4223.4403,
             4224.2413,4224.6151,4226.2994,4226.7267,4226.9876,4227.3872,4228.1580,
             4228.4166,4229.1477,4229.4545,4229.8644,4230.4267,4230.8242,4233.2873,
             4235.4636,4236.0471,4237.2198,4239.7873,4240.5947,4241.0948,4241.9725,
             4242.7197,4243.0195,4243.2611,4243.9241,4247.5988,4247.9886,4248.3908,
             4249.4761,4249.6791,4250.3146,4250.7749,4251.1846,4253.5385,4253.8674,
             4255.2373,4256.2537,4257.4963,4257.8982,4258.5204,4259.3619,4260.3330,
             4260.9845,4261.2754,4261.4920,4261.9665,4262.6124,4263.3561,4264.1064,
             4264.3399,4266.2864,4266.5271,4267.5191,4269.0639,4269.9426,4272.1689,
             4272.8745,4273.3574,4274.0246,4275.1596,4276.8073,4277.3139,4277.5282,
             4278.3232,4279.0627,4280.5680,4281.0678,4281.4145,4282.0413,4282.8976,
             4283.5184,4284.9749,4285.1823,4286.2285,4287.0808,4288.6688,4289.6551,
             4290.3937,4291.8098,4292.3057,4293.7698,4294.7193,4295.5849,4295.8155,
             4297.3066,4299.8393,4300.1008,4300.6495,4302.5261,4303.9892,4304.9565,
             4306.3668,4307.1762,4308.1220,4308.6001,4309.2392,4311.7994,4312.9975,
             4313.3084,4313.6025,4314.3195,4315.2543,4315.9484,4316.1094,4316.3726,
             4317.8400,4318.4157,4319.0973,4320.1264,4325.2741,4325.8599,4327.7145,
             4328.6879,4328.9154,4330.4120,4330.8438,4331.1995,4332.0297,4332.3395,
             4333.5612,4333.9369,4335.3379,4337.0708,4337.2774,4338.1078,4339.8716,
             4340.8954,4342.4442,4343.3817,4343.6035,4343.9515,4344.3265,4345.1680,
             4346.4367,4347.6385,4348.0640,4348.5982,4349.0722,4350.2717,4351.2719,
             4352.2049,4352.6122,4352.8201,4353.4487,4354.4824,4355.3205,4355.5259,
             4356.0445,4357.6132,4358.3200,4359.3719,4360.1672,4360.7178,4361.3071,
             4362.0662,4362.4721,4363.7945,4364.0401,4365.9301,4366.9657,4367.4180,
             4367.8316,4369.8756,4370.7532,4371.3290,4372.4900,4374.1239,4374.7851,
             4375.9542,4376.5309,4378.1768,4379.6668,4380.2863,4381.4018,4381.8601,
             4385.0566,4385.7535,4388.1014,4388.9680,4391.1105,4391.3893,4392.9740,
             4393.7590,4396.1394,4396.4783,4397.0097,4397.9149,4399.4060,4400.0968,
             4400.3863,4400.9863,4401.5810,4402.2452,4402.6001,4402.9271,4408.8828,
             4409.8991,4413.3745,4413.6324,4414.4863,4414.7670,4416.2370,4416.8447,
             4418.4157,4419.7692,4420.2571,4421.5439,4422.0479,4423.7202,4424.8373,
             4426.0011,4428.5536,4429.6348,4429.8227,4430.1890,4430.9963,4432.2523,
             4432.9628,4433.4907,4433.8380,4434.8266,4434.9564,4435.6787,4436.0484,
             4436.2844,4436.9961,4438.7464,4439.1237,4439.4614,4439.8793,4440.5737,
             4440.8659,4441.6084,4443.0863,4443.6657,4445.0338,4445.9012,4447.2296,
             4447.8346,4448.8792,4449.5206,4450.8040,4452.5655,4454.1392,4454.7738,
             4456.7080,4457.2362,4458.0015,4458.7384,4459.5477,4460.5574,4461.5278,
             4463.6660,4463.8346,4464.1386,4465.3406,4465.9751,4468.3208,4469.1891,
             4469.5253,4470.9906,4474.7594,4475.2214,4478.4032,4478.5958,4479.6376,
             4480.2661,4480.5453,4481.8107,4483.3468,4486.8973,4487.4955,4488.3123,
             4488.6803,4489.6645,4490.6715,4490.9816,4492.1007,4493.3337,4494.6907,
             4494.9616,4495.2378,4496.3155,4496.6381,4497.9144,4498.5384,4498.9401,
             4499.7037,4499.9832,4500.5769,4502.9268,4504.6306,4505.2167,4506.4730,
             4507.8339,4510.5259,4510.7332,4512.4832,4513.2231,4513.6799,4514.0573,
             4514.5503,4515.1182,4515.6274,4517.0352,4519.2592,4521.1939,4521.7407,
             4522.3230,4522.7839,4524.1289,4524.8379,4525.4456,4526.0286,4527.7086,
             4528.0208,4528.4097,4528.9740,4530.3191,4530.5523,4531.7134,4532.2577,
             4533.0770,4533.3042,4534.1198,4534.4025,4535.2546,4535.4903,4535.9798,
             4537.6426,4540.9990,4542.0513,4543.4153,4543.8692,4544.5142,4545.0519,
             4545.9156,4546.6757,4547.2498,4547.7589,4548.3261,4551.4736,4552.1536,
             4553.8376,4555.8127,4556.3872,4558.3458,4559.0571,4559.3120,4561.3477,
             4562.0803,4562.4756,4563.6609,4564.4054,4564.8354,4567.2403,4568.1426,
             4568.5235,4570.9722,4571.2195,4571.5312,4573.7038,4574.0279,4574.3058,
             4575.2338,4575.4259,4576.7405,4577.5822,4579.3495,4579.8272,4581.1732,
             4584.3671,4588.4261,4589.1222,4589.2668,4589.8978,4592.6660,4592.9675,
             4593.6437,4595.4206,4596.0967,4596.3078,4598.7627,4599.7046,4600.8717,
             4602.8862,4603.1446,4604.6529,4605.8237,4606.5021,4607.3780,4607.9347,
             4608.6201,4609.5673,4610.0809,4610.6272,4611.8597,4612.5435,4613.6044,
             4615.0240,4615.3340,4615.8838,4616.4523,4617.9065,4619.4794,4620.2408,
             4620.4476,4621.1629,4622.1209,4623.1169,4623.8911,4624.1365,4624.3140,
             4624.6776,4625.0538,4625.2747,4625.4938,4627.2980,4628.2015,4628.4409,
             4631.7617,4633.7657,4637.2328,4638.1259,4638.6849,4639.5188,4639.7036,
             4640.0462,4643.5566,4644.7072,4646.6861,4647.2509,4647.9564,4649.9777,
             4650.2343,4651.5552,4651.9895,4654.0057,4655.2126,4656.1849,4656.5518,
             4657.9012,4659.5698,4660.3941,4662.5936,4663.2026,4664.9704,4666.0050,
             4666.5158,4666.7985,4668.1716,4669.6658,4669.9842,4672.4372,4673.6609,
             4675.3760,4676.0555,4676.9721,4678.2352,4680.2377,4680.6460,4682.2319,
             4682.7322,4683.3517,4686.1946,4686.5715,4690.3343,4690.6219,4691.0516,
             4691.3452,4691.6354,4691.8842,4694.0914,4695.0381,4695.4542,4698.2248,
             4700.1402,4700.7711,4702.3161,4703.3581,4703.9898,4705.2986,4705.7606,
             4706.2511,4707.0450,4708.2940,4710.8238,4712.0056,4712.4814,4712.8408,
             4714.6715,4715.4308,4717.3905,4719.4424,4719.9795,4720.4586,4720.7811,
             4721.2765,4721.5910,4722.0886,4723.4382,4723.7840,4724.7723,4726.8683,
             4727.8491,4728.1330,4729.1282,4729.8795,4730.6593,4730.8815,4732.0532,
             4734.0458,4735.9058,4737.9172,4739.6764,4740.5292,4740.9585,4741.3039,
             4742.1174,4742.5671,4743.6870,4745.3375,4747.6176,4748.5911,4749.2002,
             4749.9713,4752.4141,4757.2196,4761.1101,4762.5242,4764.3463,4764.8646,
             4765.5950,4766.6006,4768.0578,4773.2410,4775.3130,4775.7940,4776.7792,
             4777.1915,4778.2940,4779.7286,4780.7505,4781.2902,4781.5260,4782.7612,
             4783.6389,4783.8617,4784.0396,4786.5310,4787.1480,4788.6781,4789.3868,
             4790.4373,4791.5157,4792.0819,4793.2446,4793.9052,4795.9131,4796.2429,
             4800.1724,4801.0517,4802.6714,4803.9559,4805.6063,4806.0205,4808.1337,
             4809.6140,4811.8710,4812.3755,4813.0069,4813.7204,4813.8963,4817.0206,
             4818.6477,4819.1930,4820.4649,4820.8847,4821.5878,4821.8646,4822.0140,
             4822.8548,4823.1823,4823.6058,4823.9967,4824.6484,4826.7004,4828.6596,
             4829.7973,4831.1213,4831.5975,4832.4233,4832.8025,4833.6759,4837.4121,
             4838.3546,4840.4744,4840.8492,4842.1676,4843.9413,4844.1653,4844.7553,
             4844.9251,4845.1626,4847.3263,4847.8095,4848.3625,4849.1399,4849.8617,
             4850.4397,4852.8685,4853.6087,4856.2756,4857.5385,4858.0944,4858.3327,
             4861.2167,4861.7173,4863.1724,4865.4775,4865.9100,4867.5560,4868.2742,
             4868.5277,4868.8814,4871.2890,4872.0309,4872.9169,4874.3645,4876.4949,
             4877.0015,4877.8128,4878.0094,4878.7330,4879.3497,4879.8635,4881.2046,
             4881.8531,4882.2432,4885.9843,4886.8661,4888.2612,4889.0422,4889.4903,
             4889.8553,4890.4582,4891.0379,4891.6606,4892.7592,4893.4453,4894.9551,
             4898.8044,4899.2401,4902.0545,4902.7943,4904.7516,4907.2093,4909.8431,
             4910.1576,4910.5488,4910.7929,4911.1358,4911.3787,4912.0428,4912.5293,
             4913.4789,4914.1210,4915.4155,4917.8220,4919.8157,4921.6134,4922.9442,
             4924.4223,4925.4254,4925.9501,4926.2529,4927.2988,4927.7803,4929.0860,
             4929.9850,4933.2091,4933.8521,4934.3306,4936.7746,4937.8295,4938.5047,
             4939.0605,4939.2707,4939.6422,4941.4160,4941.9174,4943.0642,4943.7076,
             4945.4587,4946.6637,4947.5752,4950.2513,4950.6264,4952.6908,4953.5320,
             4957.2963,4958.0988,4958.7240,4961.3814,4961.7260,4962.0621,4963.1881,
             4965.0795,4965.7315,4967.0628,4967.6529,4968.7556,4970.0785,4972.1597,
             4975.0917,4975.9485,4976.5938,4977.3906,4980.1859,4980.7561,4980.9513,
             4982.4875,4983.5324,4985.3725,4985.9486,4987.1470,4987.5582,4989.3086,
             4992.1257,4992.6372,4993.7488,4994.1058,4997.8189,4999.9398,5000.2463,
             5002.0972,5003.5981,5004.1279,5008.1897,5009.3344,5009.9367,5010.4174,
             5011.4774,5013.1647,5014.7539,5015.8893,5016.5356,5017.1628,5017.5084,
             5018.0535,5019.8062,5020.5458,5021.2531,5022.0051,5022.1740,5023.7084,
             5028.6556,5029.0119,5029.6295,5029.6370,5029.8916,5035.3375,5036.7284,
             5038.3018,5039.2303,5039.5265,5040.6805,5041.1224,5041.5998,5043.5145,
             5044.7195,5045.2481,5046.3527,5046.6372,5047.0434,5048.9366,5049.7960,
             5050.7842,5051.8887,5055.3473,5057.9866,5059.8611,5061.6562,5062.0371,
             5062.9325,5063.5157,5063.9988,5064.6020,5064.9454,5066.1355,5066.7773,
             5067.1379,5067.9737,5069.3384,5074.6465,5075.4659,5079.1374,5079.9079,
             5081.4462,5084.6994,5084.9935,5085.2955,5086.1774,5089.2192,5090.0513,
             5090.5455,5095.0639,5096.4848,5098.0432,5100.6211,5101.1299,5109.7331,
             5111.2781,5115.0448,5117.2923,5118.2023,5122.4995,5125.4895,5125.9502,
             5128.4897,5130.2338,5131.0720,5133.1051,5134.7460,5136.1210,5137.4733,
             5140.7736,5141.7827,5143.2673,5143.9165,5145.3083,5146.0557,5148.2115,
             5149.2069,5151.6120,5154.2430,5158.6041,5161.5396,5162.2846,5163.4584,
             5165.7728,5166.6553,5168.9225,5170.2227,5173.6715,5174.7988,5175.3248,
             5175.9115,5176.2292,5176.4037,5176.9610,5177.6228,5178.4800,5180.7209,
             5182.5269,5183.0897,5183.6044,5183.9896,5184.4534,5184.7315,5186.4132,
             5187.3374,5187.7462,5190.8720,5193.8256,5194.4576,5195.8136,5197.2360,
             5198.7999,5199.1637,5202.0087,5203.8479,5205.1522,5207.8015,5209.7246,
             5211.2305,5213.3492,5216.5966,5216.8139,5218.5271,5219.1099,5220.7068,
             5220.9263,5221.2710,5228.2246,5228.9951,5231.1597,5233.2254,5234.1069,
             5238.8137,5239.5520,5240.1968,5247.6547,5250.8727,5253.4441,5254.4648,
             5258.3602,5260.1041,5261.4720,5264.2332,5265.5514,5266.7103,5269.7927,
             5270.2660,5272.6439,5272.9270,5273.1317,5274.1188,5276.4089,5277.1460,
             5277.5002,5280.0855,5280.3441,5281.0687,5281.6285,5286.8870,5289.8974,
             5291.8163,5294.3971,5295.0886,5296.2787,5297.7431,5298.2824,5300.5234,
             5301.4043,5306.9864,5307.4655,5309.6093,5310.2665,5312.0018,5312.5288,
             5312.9045,5315.2278,5317.4944,5318.4268,5320.7699,5322.8988,5325.1438,
             5325.4317,5326.2774,5326.9756,5329.3745,5329.6979,5329.7559,5330.0804,
             5337.0154,5340.4982,5343.5812,5345.3117,5346.3778,5347.9714,5349.0054,
             5349.4608,5351.1265,5351.8365,5354.6016,5355.6365,5358.7059,5359.8269,
             5360.1501,5361.1556,5362.5751,5365.4705,5369.2819,5369.4470,5370.7096,
             5372.7027,5373.4943,5374.8219,5375.3526,5375.7688,5376.1304,5376.7809,
             5378.8356,5379.1105,5381.3709,5382.9276,5384.0356,5386.6107,5388.0507,
             5392.5726,5393.5995,5393.9719,5394.7608,5397.5160,5398.7015,5398.9219,
             5399.1745,5399.6218,5400.1455,5402.6048,5403.1993,5406.7553,5407.3439,
             5407.6535,5410.7687,5415.5226,5417.4858,5419.1285,5421.3517,5421.6628,
             5421.8359,5424.0079,5425.6783,5426.4075,5429.1048,5431.1120,5432.7650,
             5433.2919,5433.7005,5434.1514,5435.1229,5435.7236,5435.8925,5437.3877,
             5438.1374,5439.9891,5440.6015,5441.2145,5442.2477,5443.1187,5443.6893,
             5447.1536,5448.2719,5449.4788,5451.6520,5452.2187,5455.1524,5458.9676,
             5461.7358,5462.6129,5464.2055,5467.1704,5470.7591,5477.5522,5479.0751,
             5480.8907,5484.1467,5487.8414,5488.6287,5488.9254,5489.0885,5490.1194,
             5492.6435,5493.2040,5494.3306,5495.8738,5496.1370,5498.1841,5499.2554,
             5499.6478,5501.2814,5501.9439,5504.3018,5506.1128,5506.9196,5507.5385,
             5508.5585,5509.9938,5514.8731,5518.9893,5519.2353,5521.7532,5524.5824,
             5524.9570,5527.2953,5528.2272,5529.0972,5530.0763,5530.6977,5535.9710,
             5537.5563,5538.6087,5539.2618,5539.9106,5541.1466,5541.5815,5541.9357,
             5542.8900,5545.0495,5545.8058,5546.1204,5547.1339,5548.1758,5551.3719,
             5552.6229,5555.5312,5556.1487,5557.0454,5557.9214,5558.3422,5558.7020,
             5559.8912,5564.2011,5567.9986,5570.9268,5571.1913,5572.0951,5572.4649,
             5573.3535,5574.9153,5576.2045,5577.6845,5579.3583,5580.0771,5580.7547,
             5581.9668,5582.8665,5583.7621,5585.8472,5587.0263,5587.7351,5588.7506,
             5590.1127,5593.6129,5594.4604,5595.0635,5597.4756,5598.4794,5599.6548,
             5601.1216,5601.6032,5602.8630,5604.5154,5606.3861,5606.7330,5609.5735,
             5610.6809,5612.0679,5612.6174,5615.3195,5616.6908,5617.2897,5617.6417,
             5619.9755,5624.9129,5625.6782,5630.2969,5633.2950,5635.8792,5639.7461,
             5641.5575,5641.7342,5645.5249,5645.6688,5645.8947,5646.4514,5647.7069,
             5648.6863,5648.9910,5650.7043,5652.9020,5654.0234,5654.5185,5655.4897,
             5657.9255,5659.1272,5660.6594,5663.0420,5664.6210,5665.1799,5665.6283,
             5666.4185,5667.1281,5673.8358,5674.9864,5677.0529,5679.0050,5681.9001,
             5685.1921,5687.3489,5690.6934,5691.6612,5696.3904,5698.2937,5700.9176,
             5702.6511,5705.6371,5707.1033,5708.6794,5715.7244,5717.1711,5719.6227,
             5720.1828,5724.2374,5724.4638,5725.3885,5727.7103,5729.9823,5736.0297,
             5739.5196,5741.1705,5741.8290,5742.0838,5747.1931,5747.3849,5748.7412,
             5749.3883,5749.7862,5753.0265,5760.5508,5762.7941,5763.5290,5767.7785,
             5768.1812,5771.7601,5772.1143,5773.9463,5777.4006,5779.8352,5789.6451,
             5792.4304,5796.0683,5797.3194,5798.4780,5800.8297,5802.0798,5804.1412,
             5805.7017,5807.6814,5808.6569,5812.9725,5815.4219,5818.4495,5819.6027,
             5822.7933,5826.2722,5829.1097,5829.8625,5830.8272,5832.3705,5834.2633,
             5838.9502,5839.8514,5840.6404,5842.0500,5843.8069,5844.7909,5845.9188,
             5852.6806,5853.4745,5854.1207,5855.4986,5857.4491,5860.3103,5863.7184,
             5866.8119,5868.3745,5869.8507,5870.5526,5871.1828,5872.6028,5874.0383,
             5874.3513,5874.9865,5879.1265,5882.6242,5884.0331,5885.7016,5886.5314,
             5888.2622,5888.5841,5889.9531,5891.4510,5894.6979,5895.2818,5896.3631,
             5899.8441,5902.6021,5904.1592,5905.1763,5905.5707,5911.2296,5912.0853,
             5913.3638,5914.3856,5914.6709,5916.7282,5918.9442,5925.4036,5927.1258,
             5928.8130,5929.9343,5932.4274,5936.3861,5937.1619,5937.6634,5938.4570,
             5938.8252,5940.5653,5942.3585,5943.4653,5944.6472,5946.2347,5948.7991,
             5954.5852,5955.5616,5956.2588,5960.7806,5964.4723,5968.3199,5969.7373,
             5971.6008,5972.4351,5973.6649,5975.0648,5976.2005,5981.8846,5982.0960,
             5986.2665,5987.3016,5989.0447,5991.0071,5993.4941,5994.1287,5996.6297,
             5998.9987,6001.2033,6001.7339,6005.1650,6005.7242,6007.0722,6010.1606,
             6011.5337,6013.6777,6015.4218,6016.3589,6018.9934,6020.2715,6021.0357,
             6021.4113,6023.2242,6025.1500,6026.0448,6028.2713,6029.2257,6029.6500,
             6030.4450,6032.1274,6032.8726,6035.1923,6036.7622,6037.6975,6038.6805,
             6042.5898,6043.2233,6044.4327,6046.8977,6049.0510,6050.9818,6052.7229,
             6053.3808,6055.5937,6058.1813,6059.3725,6061.5361,6064.7508,6065.7799,
             6067.0803,6069.0207,6070.3412,6073.1036,6077.1057,6077.8728,6078.4212,
             6079.2227,6081.2433,6085.3748,6085.8797,6087.2622,6088.0305,6090.4482,
             6097.1945,6098.1205,6098.8031,6099.0831,6099.9886,6101.7251,6105.6351,
             6107.5337,6112.8375,6113.4657,6114.5392,6114.9234,6116.1661,6119.6998,
             6120.5564,6121.4075,6122.2144,6123.3619,6123.8331,6124.0687,6124.4805,
             6125.7396,6127.4160,6129.5458,6133.8143,6137.9260,6138.6465,6144.7581,
             6145.4411,6150.6830,6151.9929,6154.0682,6154.5162,6155.2385,6155.5810,
             6157.0878,6161.3534,6162.1700,6162.9699,6164.4796,6165.1232,6169.8221,
             6170.1740,6172.2778,6173.0964,6178.4315,6180.7050,6182.6217,6187.1350,
             6188.1251,6189.1449,6191.9053,6193.8563,6198.2227,6200.4325,6201.1002,
             6203.4925,6205.8599,6207.2201,6207.7504,6208.6869,6215.9383,6220.0112,
             6221.3192,6224.5272,6226.3697,6232.9745,6234.1962,6234.8554,6239.7121,
             6240.9536,6243.1201,6245.0409,6245.6182,6246.1736,6248.4055,6250.4858,
             6257.4237,6258.6068,6261.4181,6264.7139,6265.6036,6266.1737,6268.2001,
             6271.5444,6274.1166,6276.1646,6277.2385,6278.6453,6279.1660,6279.9776,
             6285.2780,6286.0453,6287.2554,6289.4882,6291.1915,6291.6286,6292.0601,
             6292.8909,6293.2424,6296.8722,6298.9027,6300.9165,6301.4128,6303.2507,
             6304.2423,6307.6570,6309.1600,6310.8101,6312.6222,6314.2693,6315.7750,
             6317.1824,6321.3004,6321.8199,6324.3978,6326.3669,6327.2778,6331.4137,
             6332.5118,6333.1459,6335.4060,6337.6201,6339.6684,6342.8595,6346.1209,
             6348.2321,6348.7375,6355.6300,6355.9108,6357.0229,6357.6779,6358.6211,
             6359.1345,6359.6744,6362.2518,6364.8937,6369.1394,6369.5748,6371.9436,
             6376.9305,6378.7604,6379.6733,6381.3598,6381.7588,6384.7169,6387.3957,
             6388.8134,6389.3900,6390.1289,6391.1774,6392.3685,6393.7972,6394.0498,
             6394.7289,6396.9499,6399.2065,6399.9498,6400.6961,6403.0128,6406.4462,
             6411.8991,6413.6145,6415.5374,6416.3071,6418.3703,6422.1075,6422.8969,
             6424.8127,6428.7733,6431.5550,6436.6699,6437.7613,6438.3095,6438.9160,
             6439.0715,6441.8994,6443.8598,6445.6407,6446.7712,6450.0060,6450.9557,
             6452.0591,6455.2647,6457.2824,6460.4750,6461.0596,6467.6754,6471.2131,
             6483.0825,6485.3759,6487.4811,6488.8832,6490.7372,6493.1975,6493.7777,
             6493.9694,6495.2573,6496.0292,6497.4919,6499.1061,6499.6456,6500.6575,
             6501.3608,6501.9919,6503.5111,6506.9868,6508.3605,6509.0503,6512.3639,
             6513.8462,6522.0432,6522.4994,6531.3418,6534.6056,6536.1441,6537.1775,
             6537.6139,6538.1120,6542.0498,6545.7186,6550.1877,6551.7055,6554.1603,
             6558.8756,6560.0574,6564.4445,6565.0700,6569.6322,6576.1221,6577.2146,
             6577.6552,6580.2299,6583.9060,6584.6130,6585.6945,6588.5396,6591.4845,
             6593.4628,6593.9391,6596.1141,6598.6778,6599.4824,6600.7311,6602.7620,
             6603.6209,6604.8534,6605.4165,6613.3740,6617.0580,6617.5154,6618.1668,
             6619.9458,6620.9665,6632.0837,6633.4579,6638.2207,6638.9119,6639.7403,
             6643.6976,6644.6635,6646.5407,6648.4954,6648.9586,6654.3675,6655.4892,
             6656.9386,6658.6774,6660.6761,6662.2686,6663.6966,6664.0510,6666.3588,
             6668.8163,6673.5797,6674.6969,6677.2817,6678.7067,6680.0811,6683.3673,
             6684.2929,6687.5207,6692.7262,6694.0068,6694.4967,6696.1400,6697.7123,
             6700.7499,6711.2521,6713.9700,6717.3851,6719.1995,6722.8899,6726.3091,
             6727.4583,6728.1183,6728.7595,6729.9328,6732.6491,6733.7488,6738.1800,
             6741.9212,6742.8845,6746.1355,6751.4267,6752.8335,6753.6598,6756.4528,
             6757.1092,6758.2035,6766.6117,6770.1069,6772.1745,6773.0973,6778.3123,
             6779.3242,6780.1252,6780.4131,6787.7364,6788.8403,6791.2351,6795.7982,
             6798.4868,6800.4645,6807.3190,6808.5311,6809.1000,6809.3173,6809.5094,
             6809.6773,6810.5494,6812.4915,6812.7752,6815.6123,6823.5084,6824.6774,
             6827.2488,6829.0355,6829.3472,6832.8927,6834.9246,6839.2949,6846.4790,
             6851.8842,6852.3545,6853.5234,6854.1093,6854.5128,6855.3150,6855.6891,
             6861.2688,6862.8685,6863.5350,6866.3666,6866.7634,6868.4508,6871.2891,
             6874.7530,6879.5824,6882.8114,6886.4082,6887.0881,6888.1742,6889.3032,
             6892.2502,6894.2303,6900.7632,6902.2105,6908.9905,6909.8491,6911.2264,
             6914.7127,6916.1289,6920.0394,6925.0094,6936.6528,6937.6642,6942.5380,
             6943.6105,6945.4902,6946.2134,6951.4776,6952.9666,6954.6560,6955.3150,
             6960.2500,6962.3117,6965.4307,6969.2976,6981.0831,6985.4718,6986.0298,
             6989.6553,6992.2126,6992.6966,6993.0371,6993.9873,6996.7573,6999.6238,
             7000.8036,7002.8829,7007.0961,7014.9689,7015.3172,7018.5675,7020.4841,
             7021.2832,7025.2251,7026.4616,7028.0228,7030.2514,7032.9310,7036.2831,
             7038.7202,7045.7970,7053.6196,7054.4177,7058.4895,7059.5253,7060.0417,
             7060.6538,7061.3942,7064.4515,7067.2181,7068.7358,7071.0942,7071.4800,
             7072.3939,7074.2563,7075.3336,7077.3690,7084.1690,7086.7044,7088.8232,
             7089.3395,7100.5144,7107.4778,7109.1819,7109.8602,7112.9189,7113.9879,
             7114.3984,7124.5607,7125.8200,7130.1846,7130.7236,7131.3586,7132.6100,
             7140.4617,7142.3310,7147.0416,7148.5589,7150.2844,7153.5877,7156.9416,
             7158.8387,7159.9470,7162.5569,7167.2028,7168.8952,7170.3607,7173.3725,
             7176.1860,7176.7213,7188.5320,7191.1328,7200.0454,7201.8091,7202.1910,
             7202.5166,7206.4830,7206.9804,7208.0063,7212.6896,7218.0542,7219.1515,
             7220.9809,7225.1099,7229.9380,7229.9386,7230.8623,7233.5365,7240.1848,
             7242.0919,7244.6965,7246.1277,7250.5895,7253.6760,7255.3541,7256.9860,
             7258.1770,7265.1724,7270.6636,7272.9359,7284.9033,7285.4437,7296.2659,
             7298.1434,7305.4043,7311.7159,7315.0661,7316.0050,7323.2106,7324.8073,
             7326.1491,7328.2850,7329.4916,7335.5772,7339.6035,7341.1515,7342.5769,
             7345.3996,7346.3428,7348.0530,7350.8140,7353.2930,7358.3449,7361.3472,
             7370.8255,7372.1184,7376.8772,7380.4263,7383.9805,7385.5006,7392.9801,
             7393.4377,7402.2521,7410.9683,7411.7363,7412.3368,7417.7908,7418.5499,
             7422.3118,7423.8044,7425.2942,7428.9405,7430.2533,7435.3683,7436.2970,
             7440.4933,7444.7489,7447.8488,7455.2080,7458.7541,7461.8746,7462.9910,
             7469.0566,7471.1641,7481.3545,7483.6256,7484.3267,7487.9738,7493.4276,
             7495.5641,7499.0024,7500.6556,7503.8691,7508.4787,7510.4082,7511.3498,
             7511.7902,7514.6518,7518.7824,7523.1347,7525.5079,7528.4892,7531.1437,
             7536.4365,7537.4287,7549.3138,7555.3254,7565.8515,7566.5296,7567.7417,
             7569.5115,7570.5590,7571.0326,7573.3437,7589.3151,7598.2054,7607.8230,
             7618.3443,7620.0773,7625.7053,7627.1749,7628.8818,7630.3106,7632.1305,
             7635.1060,7637.3852,7638.7805,7647.3794,7651.7451,7652.3202,7653.8284,
             7654.6998,7658.3246,7660.0234,7660.8903,7666.5679,7668.9608,7670.0575,
             7672.2552,7676.2195,7678.1267,7683.0189,7685.3075,7691.7706,7693.8016,
             7699.4034,7699.7794,7701.1081,7703.6850,7704.8169,7710.2692,7712.4050,
             7713.9378,7721.2023,7723.7611,7724.2072,7728.9510,7731.7385,7742.5628,
             7743.9552,7761.7133,7762.7320,7771.9468,7776.6727,7778.4922,7782.3165,
             7787.8002,7788.9342,7798.3579,7810.6245,7813.4761,7813.9729,7814.3230,
             7816.1534,7817.7669,7822.3877,7834.4578,7835.6209,7836.4597,7840.2981,
             7841.7911,7842.2656,7847.5394,7848.4547,7861.9100,7861.9110,7864.0218,
             7865.9698,7868.1946,7872.6316,7875.4626,7877.0756,7886.2830,7891.0750,
             7893.6191,7896.4528,7899.6229,7900.3200,7916.4420,7924.9914,7925.7478,
             7937.7335,7941.7259,7948.1764,7954.5923,7956.9740,7972.5960,7974.1592,
             7978.9731,7981.2262,7987.9731,7991.3655,7993.6808,8000.0449,8006.1567,
             8014.7857,8017.5275,8022.2014,8024.2530,8025.7271,8026.2135,8030.2004,
             8032.4313,8037.2183,8046.1169,8050.6480,8053.3085,8054.5355,8062.6304,
             8075.6518,8085.2190,8089.4852,8092.2373,8093.6238,8094.0559,8103.6931,
             8115.3110,8119.1811,8122.7234,8129.4051,8138.4753,8139.9032,8143.1380,
             8149.7021,8152.3818,8157.5434,8159.7277,8163.1210,8166.4477,8169.7865,
             8177.1788,8186.9113,8190.8849,8194.3943,8198.4420,8202.1469,8207.4785,
             8209.4461,8214.1472,8217.2264,8231.4069,8234.9381,8252.3936,8253.6156,
             8254.7422,8259.5110,8261.0143,8264.5225,8275.6266,8288.4143,8292.5272,
             8295.5497,8297.1765,8304.4244,8311.6306,8320.8554,8330.4494,8335.7067,
             8341.4770,8356.0627,8356.0693,8358.7260,8360.4915,8366.0741,8367.3936,
             8369.3404,8379.2261,8384.7240,8385.7280,8387.1053,8388.5363,8398.1780,
             8399.2578,8401.9890,8403.7965,8408.2096,8411.9172,8416.7269,8417.9982,
             8421.2254,8424.6475,8432.4935,8445.4870,8446.5116,8450.6754,8456.3468,
             8464.2367,8465.6706,8471.8260,8478.3580,8490.3065,8500.6797,8510.6240,
             8511.9091,8516.5542,8521.4422,8539.7930,8542.0871,8543.7225,8544.5955,
             8554.9440,8556.3250,8558.4464,8560.4348,8568.2088,8573.1205,8575.3305,
             8577.2768,8587.6342,8591.8387,8593.1089,8599.7553,8604.0163,8605.7762,
             8616.2219,8620.4602,8621.3225,8629.1415,8631.3565,8638.3623,8639.4416,
             8645.3087,8649.1490,8655.8760,8662.1372,8665.4855,8667.9442,8675.3983,
             8678.4083,8687.8480,8701.1209,8703.7025,8704.8601,8707.3592,8709.2341,
             8710.4142,8712.8528,8713.6547,8719.6290,8721.6595,8722.4577,8723.7177,
             8724.3761,8732.4240,8734.0234,8739.7816,8748.0309,8749.1697,8751.2062,
             8758.2434,8760.4496,8761.6862,8766.7450,8771.8602,8772.8052,8775.5733,
             8782.7156,8784.5621,8790.3761,8792.0572,8798.1720,8799.0875,8804.5894,
             8810.2540,8816.1728,8817.7431,8829.6938,8841.1822,8842.0744,8849.9165,
             8852.7916,8854.9079,8868.8334,8875.2324,8889.1939,8892.9865,8899.2971,
             8905.6581,8907.0331,8910.8566,8920.2012,8928.0925,8941.6608,8949.1227,
             8955.8467,8957.9860,8962.1468,8967.6403,8969.8667,8985.2808,8987.4079,
             8990.8935,8995.1893,8997.8763,9008.4636,9009.8832,9012.5263,9016.5903,
             9017.5912,9031.8194,9035.9195,9037.8938,9038.6911,9040.1229,9045.3532,
             9048.2501,9053.4858,9056.0813,9062.5631,9063.9600,9066.1115,9068.0230,
             9069.5817,9072.2785,9075.3945,9090.8186,9094.8289,9107.2269,9111.5314,
             9122.9674,9132.2739,9134.6919,9140.5559,9159.0346,9165.8950,9170.8220,
             9178.7795,9187.5654,9193.5928,9194.6385,9199.6841,9203.9617,9208.5811,
             9215.9197,9224.4992,9227.5119,9232.4960,9233.8574,9239.3261,9250.5782,
             9260.3254,9263.6930,9266.2070,9266.9198,9267.6895,9270.1501,9276.2732,
             9289.5624,9291.5313,9294.9740,9300.0131,9307.8956,9310.4440,9317.7296,
             9336.1509,9336.1624,9340.7053,9354.2198,9360.9879,9383.2722,9388.9308,
             9390.5852,9399.0891,9409.3488,9414.0886,9420.4985,9431.5996,9432.2827,
             9436.8127,9455.2023,9461.0279,9467.1954,9470.6819,9474.8793,9486.9256,
             9494.0442,9495.4979,9497.1891,9500.2999,9505.3930,9507.6525,9508.4513,
             9535.6566,9548.0303,9561.2452,9577.3475,9582.8133,9595.3906,9613.6911,
             9629.5693,9632.6439,9657.7863,9664.6983
        };
    }
}
