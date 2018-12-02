use std::f64;

pub const PI: f64 = 3.14159265358979323846;

fn MillimeterToMeter(m: f64)  -> f64 {
    ((m)/1000.0)
} 

fn Root(input: f64, n: f64) -> f64 {
    input.powf(1.0 / n)
}

fn SolveQuadratic(a: f64, b: f64, c: f64) -> f64 {
    let D = (b * b) - (4.0 * a * c);
    
    if D < 0.0 {
        return -1.0;
    }
    if D == 0.0
    {
        return -b / (2.0 * a);
    }
    
    let x1 = (-b + D.sqrt()) / (2.0 * a);
    let x2 = (-b - D.sqrt()) / (2.0 * a);
    
    if x1 > x2 {
        return x1;
    }
    return x2;
}

// a - коэффициент избытка топлива
// Bh - количество топлива сгораемого за час в кг
fn get_m(Bh: f64, fuel: Fuel) -> f64 {
    let Gbc = 0.55 * (fuel.C / (fuel.CO2 + fuel.CO)) + 0.0021 * fuel.C + 0.0406 * fuel.H + 0.0045 * fuel.W;
    Gbc * Bh
}
   
// a - коэффициент избытка топлива
// Bh - количество топлива сгораемого за час в кг
fn get_n(Bh: f64, fuel: Fuel) -> f64  {
    let Gbb = 0.0000445 * (fuel.C / (fuel.CO2 + fuel.CO)) + 0.0000012 * fuel.C + 0.0000044 * fuel.H + 0.0000005 * fuel.W;
    Gbb * Bh
}

// q0 - распологаемое тепло
// a - коэффициент избытка топлива
// Bh - количество топлива сгораемого за час в кг
// tb - температура наружного воздуха
fn Getq0(Bh: f64, tb: f64, fuel: Fuel, heatCoefficient: HeatCoefficient) -> f64 {
    // q0 = q0' + q0''
    // q0' = Bh * K
    // q0'' = M * tb + N * tb^2
       
    Bh * fuel.K + heatCoefficient.M * tb + heatCoefficient.N * tb * tb
}

// Q2' - химические потери тепла
// a - коэффициент избытка топлива
// b0 - химическая характеристика топлива
// Bh - количество топлива сгораемого за час в кг
fn GetQ21(Bh_fact: f64, fuel: Fuel) -> f64 {
    56.9 * fuel.C * (fuel.CO / (fuel.CO2 + fuel.CO)) * Bh_fact
}
   
// Q2'' - механические потери тепла
// q0   - распологаемое тепло
// q22  - суммарная потеря от уноса, слака и провала угля в %
   fn GetQ22(q0: f64, q22: f64) -> f64 {
       q0 * q22 / 100.0
}

// Q3 - потеря тепла с уходящими газами
// a - коэффициент избытка топлива
// Bh - количество топлива сгораемого за час в кг
// q0 - распологаемое тепло
// T3 - температура газов по выходе из трубчатой части котла
/*
internal inline f64
GetQ3(f64 Bh, f64 q0, f64 T3, Fuel fuel, HeatCoefficient heatCoefficient)
{
    auto q3 = heatCoefficient.M / fuel.K * (heatCoefficient.M * T3 + heatCoefficient.N * T3 * T3) * 100.0;
    return q0 * q3;
}
*/

fn GetQ3(T3: f64, heatCoefficient: HeatCoefficient) -> f64 {
    heatCoefficient.M * T3 + heatCoefficient.N * T3 * T3
}

// Q4 - потери на внешнее охлаждение
// phi - коэффициент качества изоляции
// H0 - наружная поверхность котла м^2
// V - скорость в км/ч
// tb - температура наружного воздуха
/*
fn GetQ4(phi: f64, H0: f64, V: f64, tb: f64, tk: f64) -> f64 {
    // Q4 = 2.2 * phi * H0 * (tk - tb)^(4/3)
    return phi * H0 * (2.2 + 0.21 * pow(V, 0.7)) * pow(tk - tb, 4.0 / 3.0);
}
*/

// Q4 - потери на внешнее охлаждение
// k  - потеря в % от q0
fn GetQ4(q0: f64, k: f64) -> f64 {
    return q0 * k / 100.0;
}

// Q5 - потеря на служебные нужды
fn GetQ5(q0: f64) -> f64 {
    // Q5 = 2-5%
    q0 * 0.035
}

// k - коэффициент теплопередачи
fn GetK1(d: PipeDiameter, omega: f64) -> f64 {
    let r = d.GetHydraulicRadius();
    6.0 + 2.45* omega.powf(0.7) * (0.0115/r).powf(0.214)
}

// k - коэффициент теплопередачи
fn GetK2(heatCoefficient: HeatCoefficient,  Hd: f64, tk: f64, T2: f64, T3: f64) -> f64 {
    let a= heatCoefficient.M + 2.0*heatCoefficient.N*tk;
    let b = ((T2 - tk)/(T3-tk)).ln();
    let c = 2.0*heatCoefficient.N*(T2-T3);
    (a * b + c) / Hd
}

//v - удельный объем протекающих по трубкам продуктов сгорания
fn GetV(Tabs: f64) -> f64 {
    (29.27*Tabs)/10330.0
}

// omega - скорость протекания газов по дымогарным трубам
fn GetOmega(fuel: Fuel, d: PipeDiameter, Tabs: f64, L0: f64, Bh: f64) -> f64 {
   let Od = d.GetPipeSquare();
   let v = GetV(Tabs);
   let top = (L0 * fuel.alpha + 1.0) * Bh * v;
   let bottom = 3600.0 * Od;
   top / bottom
}

// Дымогарные трубы

fn GetKdHd( b: f64,  M: f64,  N: f64,  T2: f64,  Tk: f64,  tk: f64,  T3d: f64) -> f64 {
    (1.0 - b) * ((M + 2.0 * N * tk) * ((T2 - Tk) / (T3d - tk)).ln() + 2.0 * N * (T2 - T3d))
}

// Qt  - тепло проходящее в котел через топочную
fn GetQt(q0: f64, Q21: f64, Q22: f64, T2: f64, heatCoefficient: HeatCoefficient) -> f64 {
    return q0 - Q21 - Q22 - (heatCoefficient.M * T2 + heatCoefficient.N * T2 * T2);
}

// T1   - Действительная температура горения
// q0   - распологаемое тепло
// Q2'  - химические потери тепла
// Q2'' - механические потери тепла
// a    - коэффициент избытка топлива
// Bh   - количество топлива сгораемого за час в кг
fn GetT1(q0: f64, Q21: f64, Q22: f64, Bh: f64, fuel: Fuel, heatCoefficient: HeatCoefficient) -> f64 {
    //N * T^2 + M * T - Q / 0.84 = 0
    let Q = (q0 - (Q21 + Q22)) / 0.84;
    return SolveQuadratic(heatCoefficient.N, heatCoefficient.M, -Q);
}

// T2 - температура при входе газов в отверстия огневой решетки (43)
// Bh - количество топлива сгораемого за час в кг
// Ht - поверхность нагрева огневой коробки
fn GetT2(Bh: f64, Ht: f64, fuel: Fuel) -> f64 {
    const A: f64 = 1350.0;

    let top = (Bh * fuel.K) / Ht + 4400.0;
    let bottom = (Bh * fuel.K) / Ht + 223000.0;

    A * Root(top / bottom, 1.6)
}

// T3 - температура газов по выходе из трубчатой части котла
// d  - диаметр дымогарной трубы
// L  - длина трубчатой части котла
// Bh - количество топлива сгораемого за час в кг
// Hk - поверхность нагрева котла
// Hi - поверхность нагрева пароперегревателя
fn GetT3(d: PipeDiameter ,  n: u64,  L: f64, Bh: f64, Hk: f64, fuel: Fuel) -> f64 {
    let Hi = d.GetH(n, L);
    let H = Hk + Hi;               // полная поверхность нагрева
    let r = d.GetHydraulicRadius(); // гидравлический радиус
    let A = (r / 0.0105).powf(0.15) * (710.0 + 332000.0 / (L / r + 105.0));
    let top = (Bh * fuel.K) / H + 4400.0;
    let bottom = (Bh * fuel.K) / H + 223000.0;

    A * Root(top / bottom, 1.6)
}

// L0 - теоретический расход воздуха для сжигания 1 кг топлива
fn GetL0(fuel: Fuel) -> f64 {
    1.0 / 23.6 * (8.0 / 3.0 * fuel.C + 8.0 * fuel.H + fuel.S - fuel.O)
}

// Lv - действительный расход воздуха для сжигания 1 кг топлива
fn GetLv(L0: f64, fuel: Fuel) -> f64 {
    L0 * fuel.alpha
}

// Ki -
fn GetKi(beta: f64, Bh: f64, rz2: f64, omegaZ2: f64, fuel: Fuel) -> f64 {
    let x = ((beta * Bh * fuel.u * fuel.K) / (10000.0 * omegaZ2)).powf(0.7);
    let y = (0.0115 / rz2).powf(0.214);
    6.0 + 0.24 * x * y
}

// Tabs - средняя абсолютная температура газов в дымогарных трубах
fn GetTabs(T2: f64, T3: f64) -> f64 {
    (T2 + T3) / 2.0 + 273.0
}

#[derive(Clone, Copy)]
pub struct PipeDiameter
{
    pub DOut: f64, // внешний диаметр в метрах
    pub DIn : f64, // внутренний диаметр в метрах
}

#[derive(Clone, Copy)]
pub struct FireChamber
{
    pub TopLengh: f64,
    pub TopWidth: f64,
    pub BottomLength: f64,
    pub BottomWidth: f64,
    pub FrontHeight: f64,
    pub RearHeight: f64,
    
    pub U: f64,
    pub R: f64,
    pub Ht: f64,
}

#[derive(Clone, Copy)]
pub struct Boiler
{
    pub Hdg: f64,           // поверхность нагрева дымогарных труб (газовая)
    pub Hz1: f64,           // поверхность нагрева жаровых труб (газовая)
    pub Hz2: f64,           // поверхность нагрева жаровых труб (газовая)
    pub Hk: f64,            // полная испаряющая поверхность нагрева котла (водяная)
    pub Hi: f64,            // поверхность нагрева пароперегревателя (газовая)
    pub l_d: f64,            // Длина труб дымогарных
    pub Lz1: f64,           // Длина труб жаровых до перегревателя
    pub Lz2: f64,           // Длина труб жаровых в области перегревателя
    
    pub Nd: u16,            // число дымогарных труб
    pub Nz: u16,            // число жаровых труб
    
    pub Dd: PipeDiameter,   // диаметр труб дымогарных
    pub Dz1: PipeDiameter,  // диаметр труб жаровых до перегревателя
    pub Dz2: PipeDiameter,  // диаметр труб жаровых в области перегревателя
    pub Di: PipeDiameter,   // диаметр труб перегревательных
    
}

#[derive(Clone, Copy)]
pub struct Fuel
{
    pub C: f64,
    pub H: f64,
    pub S: f64,
    pub O: f64,
    pub N: f64,
    pub W: f64,
    pub A: f64,
    pub CO: f64,
    pub CO2: f64,
    pub O2: f64,
    pub N2: f64,
    pub K: f64,      // теплопроизводительность 1 кг топлива (низшая)
    pub beta0: f64,     // химическая характеристика топлива
    pub alpha: f64,  // коэффициент избытка топлива
    pub u: f64,
}

#[derive(Clone, Copy)]
pub struct HeatCoefficient
{
    pub M: f64,
    pub N: f64,
}

impl PipeDiameter {

    fn new(dOut: f64, dIn: f64) -> Self {
        PipeDiameter { DOut: dOut, DIn: dIn }

    }
    // O - площадь живого сечения трубы
    fn GetPipeSquare(self) -> f64
    {
       let dSq = self.DIn * self.DIn;
       (PI * dSq) / 4.0
    }

    fn GetHydraulicRadius(self) -> f64 {
       self.DIn / 4.0
    }
}

impl HeatCoefficient {


    pub fn new(fuel: Fuel, Bh: f64) -> Self {
        let m = get_m(Bh, fuel);
        let n = get_n(Bh, fuel);
       HeatCoefficient {
           M: m,
           N: n,
       }
   }
}

impl FireChamber {

    fn new(topL: f64, botL: f64, 
           topW: f64, botW: f64,
           frontH: f64, rearH: f64, 
           u: f64) -> Self {
        FireChamber {
            TopLengh: topL,
            TopWidth: topW,
            BottomLength: botL,
            BottomWidth: botW,
            FrontHeight: frontH,
            RearHeight: rearH,
    
            U: u,
            R: botL * botW,
            Ht: topL * topW +
                rearH * topW +
                frontH * topW +
                frontH * topL +
                rearH * topL,
        }
    }

    // Bh - сжигаемое топливо в час (кг)
    // R  - площадь колосниковой решетки
    // U  - напряжение колосниковой решетки
    pub fn get_Bh(self) -> f64 {
       self.R * self.U
    }
}

impl PipeDiameter {
    // H - поверхность нагрева труб 
    // n - число труб
    // l - длина трубы
    pub fn GetH(self, n: u64, l: f64) -> f64 {
       PI * self.DOut * n  as f64 * l
    }
}

impl Fuel {

    pub fn new() -> Self {
        Fuel {
            C: 80.0,
            H: 3.1,
            S: 1.9,
            O: 3.1,
            N: 1.1,
            W: 3.1,
            A: 7.6,
            K: 7203.0,
            
            alpha: 1.35,
            
            O2: 6.2,
            N2: 79.6,
            CO: 1.5,
            CO2: 12.7,

            // beta0 - химическая характеристика топлива
            // b0 = 2.37 * (H - O/8) / C
            beta0: 2.37 * ((3.1 - 3.1 / 8.0) / 80.0),

            u: 0.0,
        }
    }

    pub fn get_Bh_fact(mut self, fire_chamber: FireChamber, Q22: f64) -> f64 {
       self.u = (100.0 - Q22) / 100.0;
       self.u * fire_chamber.R * fire_chamber.U
    }
}

fn main() {
    println!("Hello, worl_d!");


    const U: f64 = 550.0;     // напряжение колосниковой решетки
    
    const t: f64 = 0.0;     // температура окружающего воздуха °C
    const v: f64 = 0.0;       // скорость км/ч
    
    const q22: f64 = 36.0;    // суммарная потеря от уноса, шлака и провала угля в %
    const tk: f64 = 190.0;
    
    let fire_chamber = FireChamber::new(
        MillimeterToMeter(2222.0), MillimeterToMeter(2278.0), 
        MillimeterToMeter(1333.0), MillimeterToMeter(1028.0),
        MillimeterToMeter(1815.0), MillimeterToMeter(1605.0),
        550.0);
    

    
   
    /*
    fireChamber.TopLengh     = MillimeterToMeter(2649);
    fireChamber.BottomLength = MillimeterToMeter(2744);
    
    fireChamber.TopWidth     = MillimeterToMeter(1376);
    fireChamber.BottomWidth  = MillimeterToMeter(1016);
    
    fireChamber.FrontHeight  = MillimeterToMeter(1800);
    fireChamber.RearHeight   = MillimeterToMeter(1606);
    */
    
    //const f64 R = 2.8;       // площадь колосниковой решетки
    //const f64 Ht = 15.6;     // поверхность нагрева огневой коробки
    
    // дымогарная труба
    let dd = PipeDiameter::new(MillimeterToMeter(51.0), MillimeterToMeter(46.0));
    
    // длина труб дымогарных
    let l_d = MillimeterToMeter(4550.0); 

    // число дымогарных труб
    let nd : u64 = 210;                 
    
    let Hd = dd.GetH(nd, l_d);
    
    let fuel = Fuel::new();
        
    // сжигаемое топливо в час (кг)
    let Bh = fire_chamber.get_Bh();                    

    // фактически сжигаемое топливо в час (кг) 
    let Bh_fact = fuel.get_Bh_fact(fire_chamber, q22);   
    
    let heatCoefficient = HeatCoefficient::new(fuel, Bh_fact);
    
    let Q0 = Getq0(Bh, t, fuel,heatCoefficient);
    let Q21 = GetQ21(Bh_fact, fuel);
    let Q22 = GetQ22(Q0, q22);
    
    let T1 = GetT1(Q0, Q21, Q22, Bh_fact, fuel, heatCoefficient);
    let T2 = GetT2(Bh_fact, fire_chamber.Ht, fuel);
    
    //auto Q4 = GetQ4(.4, 51.9, v, t);
    let Q4 = GetQ4(q0, 1.0);                               // 1% потерь от q0
    
    let L0 = GetL0(fuel); // теоретический расход воздуха для сжигания 1 кг топлива
    
    println!("Коэффициент избытка топлива alpha\t\t {:.2}", fuel.alpha);
    println!("Химическая характеристика топлива b0 \t\t {:.2}", fuel.beta0);
    
    /*
    println!("Коэффициент уравнения тепла M \t\t\t%.2lf\n", heatCoefficient.M);
    println!("Коэффициент уравнения тепла N \t\t\t%.5lf\n", heatCoefficient.N);
    */
    
    println!("Теоретический расход воздуха L0\t\t\t{:.2} кг", L0);
    
    println!("Сжигаемое топливо в час  Bh \t\t\t{:.2} кг", Bh);
    println!("Сжигаемое топливо в час по факту \t\t{:.2} кг", Bh_fact);
    
    println!("Действительная температура горения T1\t\t{:.2} °C", T1);
    println!("Температура при входе газов в отверстия T2\t{:.2} °C", T2);
    
    /*
    const f64 betaN = 0.000001; // слой сажи в метрах
    const f64 betaS = 0.000001; // слой накипи в метрах
    
    auto TWaterSide = GetTOfWallOnWaterSide(1300.0, 190.0, betaN, betaS);
    auto TGasSide = GetTOfWallOnGasSide(1300.0, 190.0, betaN, betaS);
    
    println!("Температура стенки со стороны газов\t\t%.2lf °C\n", TGasSide);
    println!("Температура стенки со стороны воды\t\t%.2lf °C\n", TWaterSide);
    */
    
    let t_3 = GetT3(dd, nd, l_d, Bh_fact, fire_chamber.Ht, fuel);
    let Tabs = GetTabs(T2, t_3);
    
    let q_3 = GetQ3(t_3, heatCoefficient);
    
    // Qt  - тепло проходящее в котел через топочную
    let q_t = GetQt(q0, Q21, Q22, T2, heatCoefficient);
    
    println!("Температура газов на выходе котла T3\t\t{:.2} °C", t_3);
    println!("Ср. абс. температура газов в д.трубах Tabs\t{:.2} °C", Tabs);
    println!("Располагаемое тепло q0 \t\t\t\t{:.2} кал/час \t\t{:.1}%", q0, q0/q0*100.0);
    println!("Химические потери тепла Q2'\t\t\t{:.2} кал/час \t\t{:.1}%", Q21, Q21/q0*100.0);
    println!("Механические потери тепла Q2''\t\t\t{:.2} кал/час \t\t{:.1}%", Q22, Q22/q0*100.0);
    println!("Потеря тепла с уходящими газами Q3\t\t{:.2} кал/час \t\t{:.1}%", q_3, q_3/q0*100.0);
    println!("Потери на внешнее охлаждение Q4\t\t\t{:.2} кал/час \t\t{:.1}%", Q4, Q4/q0*100.0);
    
    println!("Тепло проходящее в котел через топочную Qt\t{:.2} кал/час \t\t{:.1}%", q_t, q_t/q0*100.0);
    
    let Hi = fire_chamber.Ht + Hd;
    
    println!("Площадь колосниковой решетки R\t\t\t{:.2}", fire_chamber.R);
    println!("Поверхность нагрева огневой коробки Ht\t\t{:.2} м2", fire_chamber.Ht);
    println!("Поверхность нагрева дымогарных труб Hd\t\t{:.2} м2", Hd);
    println!("Полная испаряющая поверхность Hi\t\t{:.2} м2\n", Hi);
    
    
    let omega = GetOmega(fuel, dd, Tabs, L0, Bh_fact);
    let k1 = GetK1(dd, omega);
    let k2 = GetK2(heatCoefficient,  Hd, tk, t_3, t_3);
    
    let rd = dd.GetHydraulicRadius();
    let sqr_d = dd.GetPipeSquare();
    
    //auto bb = GetBeta(l_d, rd, sqr_d, 0, 0, 0, 0, 0, 0);
    
    println!("k1 \t\t\t\t\t\t{:.2}", k1);
}
