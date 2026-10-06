# Samanyolu Simülasyonu — Yöntem ve Sonuçlar

Bu belge N-cisim koşularının yöntemini, sayısal doğrulamasını ve bulgularını sunar. Model tanımı ve
kabul kriterleri için bkz. [PROJE_TANIMI.md](../PROJE_TANIMI.md).

**Özet.** Öz-yerçekimli, 40 000 parçacıklı bir yıldız diski sabit bir karanlık madde halesi içinde
2000 Myr boyunca evriltilmiştir. Hale kütlesi yarıya indirildiğinde çubuk daha erken, daha güçlü ve
daha uzun oluşmakta, 2000 Myr boyunca kalıcı olmaktadır. Bu çubuğun eşdönme yarıçapı / uzunluk oranı
**ℛ = 1.11**'dir (faz-sabit uzunlukla; yarı-maksimum uzunlukla 1.06); bu değer gözlenen çubukların
**1.0–1.4** aralığındadır, yani çubuk "hızlı" sınıfındadır. Toomre Q = 0.8 ile başlatılan disk
parçalanmakta ve kendi kendini ısıtarak kararlı bölgeye geçmektedir. Normal haleli modelde çubuğun
1300–1800 Myr arasındaki zayıflaması test edilen dört hipotezle açıklanamamıştır (§5).

---

## 1. Yöntem

### 1.1 Koşular

Tüm koşularda θ = 0.6 kuadrupol momentli Barnes–Hut ağacı, eps = 0.08 kpc, dt = 0.2 Myr kullanılmıştır.
Snapshot aralığı 10 Myr'dir (2000 Myr'lik koşularda 201 snapshot).

| Koşu | N | Q_target | M_vir | Süre | Öz-yerçekimi | Hesap süresi |
|---|---|---|---|---|---|---|
| `ana` | 40 000 | 1.3 | 1.0e12 M☉ | 2000 Myr | ağaç | 1488 s |
| `soguk` | 40 000 | 0.8 | 1.0e12 M☉ | 2000 Myr | ağaç | 1613 s |
| `hafif_hale` | 40 000 | 1.3 | 5.0e11 M☉ | 2000 Myr | ağaç | 1508 s |
| `dogrulama_dogrudan` | 5 000 | 1.3 | 1.0e12 M☉ | 600 Myr | doğrudan N² | 55 s |
| `dogrulama_agac` | 5 000 | 1.3 | 1.0e12 M☉ | 600 Myr | ağaç (aynı seed) | 28 s |
| `dogrulama_dogrudan_seed7` | 5 000 | 1.3 | 1.0e12 M☉ | 600 Myr | doğrudan N², farklı seed (gürültü tabanı) | 54 s |

Disk kütlesi tüm koşularda M_d = 4.82e10 M☉'dir. `hafif_hale` koşusunda disk kütlesi yeniden
çözülmemiştir; bu modelde v_c(8.2) 229 km/s'den düşüktür.

### 1.2 Çubuk ölçümü

- A2/A0, 0.5 kpc'lik kuşaklarda hesaplanmıştır.
- **Çubuk gücü** = A2/A0(R)'nin R < 5 kpc'deki maksimumu. Yalnızca en az 200 parçacıklı kuşaklar
  hesaba katılmıştır.
- **Çubuk uzunluğu** = tepe kuşağından dışarı doğru m=2 fazının tepe fazından ±10° içinde kaldığı son
  kuşağın dış kenarı.
- **Çubuk var** = güç ≥ 0.2 ve uzunluk ≥ 2 kpc. Uzunluk koşulu olmadan, merkezdeki sarılmış sarmallar
  50–120 Myr'de sahte bir "çubuk oluşumu" vermektedir.
- **Oluşum zamanı** = bu koşulun 100 Myr kesintisiz sağlandığı ilk an.
- Gürültü tabanı: t = 0'da güç 0.034 (N=40k), 0.06–0.12 (N=5k).
- Tüm yarıçaplardaki maksimum sarmal kolları da yakalar; karşılaştırma için grafiklerde ince çizgiyle
  gösterilmiştir. R < 5 kpc toplu A2 ölçüsü ayrıca `a2_toplu_R5.png`'de verilmiştir.

### 1.3 Toomre Q ölçümü

Σ ve σ_R parçacıklardan, κ ise parçacıkların ve sabit bileşenlerin halkalar üzerinde ölçülen
v_c'sinden hesaplanmıştır. Disk zamanla yayıldığı için başlangıçtaki analitik κ kullanılmamıştır.
t = 0'da bu yöntem analitik κ ile R > 4 kpc'de %1 içinde uyuşmakta, içeride %4–6 düşük kalmaktadır
(kalın diskin gerçek alanı). Çubuk bölgesinde σ_R akış hareketlerini de içerir; orada Q yalnızca
kabaca anlamlıdır.

### 1.4 Parçalanma ölçütü

**Topak kesri**: R > 2 kpc'de, 0.25 kpc'lik piksellerde, aynı yarıçaptaki azimut ortalamasının
3 katından yoğun piksellerdeki kütlenin toplam kütleye oranı. Sentetik testte pürüzsüz bir diskte
0.004, topaklı bir diskte 0.27 vermiştir. N = 5000'de Poisson gürültüsü nedeniyle anlamlı değildir
(t = 0'da 0.30); yalnızca N = 40000 koşularında kullanılmıştır.

### 1.5 Desen hızı ve eşdönme

- **Çubuk fazı:** m=2 fazı, faz-sabit çubuk bölgesinden (0.5 kpc < R < L_çubuk) ölçülmüştür.
- **Desen hızı:** Ω_p, açılmış fazın ±50 Myr penceresindeki doğrusal eğimidir. m=2 fazı 180°'de bir
  tekrarladığı için açma, adım başına < 90° dönmeyi varsayar. t ≥ 500 Myr'de en büyük adım 23–41°'dir,
  yani açma güvenlidir. Daha erken dönemde 60–89°'lik adımlar bulunduğundan Ω_p yalnızca t ≥ 500 Myr
  için verilmiştir.
- **Eşdönme yarıçapı:** R_CR, o snapshot'ın gerçek alanından ölçülen v_c ile çözülmüştür:
  v_c(R_CR) = Ω_p · R_CR.
- **İkinci uzunluk tanımı:** A2/A0(R)'nin tepeden dışarı doğru tepe değerinin yarısına düştüğü yarıçap.

### 1.6 Merkezi kütle

R < 1 kpc içindeki **disk parçacıklarının** kütlesidir. Sabit şişkinlik ve kara delik değişmediği için
dahil edilmemiştir.

---

## 2. Doğrulama

### 2.1 Kabul kriterleri

| # | Ölçüt | Beklenen | Ölçülen |
|---|---|---|---|
| T1 | v_c(8.2 kpc) | 229 ± 5 km/s | **228.7 km/s** |
| T2 | v_c, 5–20 kpc | 205–240 km/s, düz | 212.1–232.1 km/s |
| T3 | 2πR/v_c | 220 ± 10 Myr | **220.3 Myr** |
| T4 | L_z korunumu | < 1e-6 | 1.3e-14 |
| T5 | Enerji, 1000 adım | < 1e-3 | 2.2e-6 |
| T6 | Hale kapalı, v_c(20 kpc) | < 165 km/s | 122.8 km/s |
| T7 | M_bh × 100, Δv_c(8.2) | < 1 km/s | 0.49 km/s |
| T8a | M_dyn(8.2 kpc) | 1.00e11 ± 0.05e11 M☉ | 9.970e10 M☉ |
| T8b | M_true(8.2 kpc) | 0.90e11 ± 0.05e11 M☉ | 8.970e10 M☉ (M_dyn/M_true = 1.1115) |
| T9 | Parçacıklardan v_c(8.2) | 229 ± 5 km/s | 227.1 km/s |

Asimetrik sürüklenme (8.2 ± 1 kpc): ölçülen v_c − ⟨v_φ⟩ = 8.89 ± 0.39 km/s, Jeans denkleminin
öngördüğü 8.90 km/s (oran 0.998).

### 2.2 Kuvvet doğruluğu

Ağaç kuvvetleri rastgele 500 parçacıkta doğrudan toplamla karşılaştırılmıştır (N = 20 000, başlangıç
diski). Değerler bağıl kuvvet hatasıdır.

| θ | Monopol RMS | Monopol en kötü | Kuadrupol RMS | Kuadrupol en kötü | Kazanç |
|---|---|---|---|---|---|
| 0.4 | 0.43 % | 0.73 % | 0.018 % | 0.08 % | 23× |
| 0.5 | 0.70 % | 1.21 % | 0.046 % | 0.18 % | 15× |
| **0.6** | 1.01 % | 1.78 % | **0.104 %** | 1.18 % | 9.7× |
| 0.7 | 1.37 % | 2.95 % | 0.166 % | 1.28 % | 8.3× |

Hedef RMS < %0.3'tür; seçilen ayar θ = 0.6 + kuadrupoldür. 600 Myr'deki gelişmiş diskte kazanç 8–22
kattır; adım süresi monopole göre %20 artar (42 → 50 ms).

### 2.3 Enerji korunumu

**Hata kaynağının ayrıştırılması** (600 Myr, monopol ağaç):

| Koşu | Maks \|ΔE/E\| | Davranış |
|---|---|---|
| Öz-yerçekimi kapalı, N = 20k | 1.8e-6 | ilk ~50 Myr'den sonra düz |
| Doğrudan N², N = 5k | 3.5e-5 | sınırlı, işaret değiştiriyor |
| Ağaç θ = 0.6, N = 5k, aynı seed | 6.0e-4 | tek yönde artıyor |
| Ağaç θ = 0.4, eps = 0.08 | 2.8e-4 | tek yönde artıyor |
| Ağaç θ = 0.6, eps = 0.12 | 4.0e-4 | tek yönde artıyor |
| Ağaç θ = 0.4, eps = 0.12 | 2.5e-4 | tek yönde artıyor |
| Doğrudan N², eps = 0.12 | 1.9e-5 | sınırlı |

Enerji muhasebesi ve zaman adımı temizdir; sapma monopol ağaç kuvvetinin korunumlu olmamasından gelmektedir.

**Monopol ve kuadrupol** (600 Myr, θ = 0.6, eps = 0.08):

| Koşu | Maks \|ΔE/E\| | 2000 Myr tahmini |
|---|---|---|
| monopol, N = 20k | 4.0e-4 | 1.3e-3 |
| kuadrupol, N = 20k | 2.0e-5 | ~1e-5 |
| monopol, N = 40k | 3.1e-4 | 1.06e-3 |
| kuadrupol, N = 40k | 1.6e-5 | ~6e-5 |

**Üretim koşuları** (2000 Myr, kuadrupol):

| Koşu | Maks \|ΔE/E\| | Maks \|r_com\| | Maks \|v_com\| |
|---|---|---|---|
| ana | **3.0e-5** | 0.15 kpc | 4.1 km/s |
| soguk | **8.1e-5** | 0.27 kpc | 12.7 km/s |
| hafif_hale | **5.8e-5** | 0.24 kpc | 5.3 km/s |

- T5 sınırı (1e-3) uzun koşularda 12–33 kat altında kalmıştır.
- **Kütle merkezi sürekli kaymamakta, salınmaktadır.** Soğuk koşudaki 0.27 kpc'lik sapma parçalanma
  evresinde (200–300 Myr) oluşmakta, 500 Myr'den sonra ~0.05 kpc'ye geri dönmektedir.

### 2.4 Yakınsama çalışması

600 Myr'lik koşular (monopol ağaç, kuadrupolün eklenmesinden önce). Değerler 400–600 Myr ortalamasıdır;
A2/A0 burada R < 5 kpc toplu ölçüdür.

| Koşu | R_d [kpc] | Kalınlık, medyan \|z\| [kpc] | A2/A0 |
|---|---|---|---|
| A: θ = 0.6, dt = 0.2, N = 20k | 2.89 ± 0.03 | 0.1925 | 0.09 ± 0.03 |
| B: θ = 0.5 | 2.84 ± 0.06 | 0.1930 | 0.11 ± 0.02 |
| C: dt = 0.1 | 2.82 ± 0.05 | 0.1929 | 0.08 ± 0.02 |
| D: N = 40k | 2.86 ± 0.06 | **0.1812** | 0.11 ± 0.04 |
| A′: A ile aynı, seed 7 (gürültü tabanı) | 2.79 ± 0.11 | 0.1910 | 0.10 ± 0.05 |

- R_d ve A2/A0 dört ayarda da seed gürültüsü içinde aynıdır.
- θ ve dt kalınlığı değiştirmemekte, N değiştirmektedir: fark 0.011 kpc, seed farkı 0.0015 kpc. Bu,
  iki-cisim ısınmasının bir işaretidir. Üretim koşuları bu nedenle N = 40000 ile yapılmıştır.
- Kuadrupol eklendikten sonra 600 Myr'deki kalınlık değişmemiştir (N = 20k: 0.203 → 0.200;
  N = 40k: 0.188 → 0.186); kalınlaşma ağaç hatasından değil, parçacık sayısından kaynaklanmaktadır.

### 2.5 Ağaç ve doğrudan toplama

| Karşılaştırma | Ortalama \|ΔR_d\| | Ortalama \|Δ(A2/A0)\| (toplu) |
|---|---|---|
| ağaç+kuadrupol − doğrudan (aynı seed) | **0.065 kpc** | **0.037** |
| doğrudan seed 42 − doğrudan seed 7 (gürültü tabanı) | 0.211 kpc | 0.041 |

Yaklaşık yöntemle tam yöntem arasındaki fark gürültü tabanından küçüktür; iki yöntem istatistiksel
olarak aynı sonucu vermektedir. Son %25'teki ortalama çubuk gücü 0.25 (ağaç), 0.37 (doğrudan), 0.29
(seed 7); uzunluk medyanları 2.5–3.0 kpc ve üçünde aynıdır. Güçteki fark seed farkıyla aynı
mertebededir.

N = 5000'de oluşum zamanı kaotiktir: aynı seed'le doğrudan toplamada 160 Myr, ağaçta 390 Myr; farklı
seed'le doğrudan toplamada 180 Myr. Bu düşük N'de oluşum zamanı yöntemden çok gürültüye bağlıdır.

---

## 3. Sonuçlar

### 3.1 Çubuk: güç, uzunluk, oluşum zamanı

Grafikler: `a2_karsilastirma.png` (üstte güç, altta uzunluk), `cubuk_faz_profilleri.png`
(A2/A0(R) ve faz(R)).

| Koşu | Oluşum | Çubuk var (tüm / son %25) | Tepe gücü | Son %25 ort. güç | Uzunluk medyanı (son %25) |
|---|---|---|---|---|---|
| ana | 390 Myr | %54 / **%14** | 0.38 @ 480 Myr | 0.15 | 3.5 (3.0) kpc |
| soguk | **140 Myr** | %98 / %100 | 0.50 @ 250 Myr | 0.24 | 2.5 (2.5) kpc |
| hafif_hale | **150 Myr** | %91 / **%100** | **0.50** @ 480 Myr | **0.35** | **5.0 (6.0) kpc** |

hafif_hale'de çubuk ana koşudan **240 Myr daha erken** oluşmakta, **daha güçlü** ve **daha uzundur**:
uzunluk medyanı 5.0 kpc, ana koşuda 3.5 kpc (~1.4 kat). Bu çubuk 2000 Myr boyunca kalıcıdır. Ana
koşuda çubuk oluşmakta, ancak son çeyrekte büyük ölçüde sönmektedir.

- **Uzunluk ölçüsü anlık olarak gürültülüdür** (2–9 kpc arası): fazı çubukla hizalanan bir sarmal kol
  çubuğu "uzatabilir". Bu nedenle medyanlar kullanılmıştır.
- **hafif_hale'de çubuk uzunluğu 5 kpc'de kesilmemektedir.** 5 kpc sınırı yalnızca tepe aramasına
  uygulanır; faz yürüyüşü 10 kpc'ye kadar gider. Çubuk olan anların %46'sında L > 5 kpc ölçülmüştür.
  Arama 10 kpc'ye açıldığında tepe yalnızca 12 snapshot'ta 5 kpc'nin dışına kaymaktadır; hepsi
  t ≤ 310 Myr'deki sarmal evresindedir. Faz çoğu anda 5–5.5 kpc'den sonra düzgünce kaymaktadır
  (çubuktan sarmala geçiş). 1700–1800 Myr'de faz 7 kpc'ye kadar sabittir, ancak oradaki A2 düşüktür
  (0.1–0.2); bu muhtemelen hizalanmış bir sarmal koldur.

### 3.2 Soğuk disk (Q = 0.8): parçalanma ve öz-düzenleme

Grafikler: `soguk/Q_profili.png`, `soguk/parcalanma.png`.

Disk ilk 200 Myr'de parçalanmaktadır: topak kesri 150–200 Myr'de 0.20'ye çıkmakta, ana koşuda aynı
anda 0.05–0.06'da kalmaktadır.

**Öz-düzenleme.** Kararsız disk kendi kendini ısıtmakta ve Q kararlı bölgeye yükselmektedir.

| t | Ortalama Q (1–12 kpc) | Aralık |
|---|---|---|
| 0 Myr | **0.80** | 0.76–0.84 |
| 200 Myr | 1.80 | 1.50–2.24 |
| 800 Myr | **2.32** | 2.05–2.92 |

800 Myr'de **her yarıçapta Q ≥ 2.0**'dır. Topaklar 500 Myr'de kaybolmuş, disk yayılmıştır
(R_d 2.6 → 4.3 kpc). Çarpışmasız yıldız diski kararsızlığı ısınarak doyurmaktadır; gaz diskinden
farklı olarak kalıcı topak oluşmamaktadır.

### 3.3 Desen hızı ve eşdönme

Grafik: `desen_hizi.png`.

| Koşu | Ω_p, 500 Myr | Ω_p, 1300 Myr | Ω_p, 1800 Myr | 1300–1800 eğilimi |
|---|---|---|---|---|
| ana | 48.6 | 41.4 | 44.1 km/s/kpc | **+0.56 /100 Myr** (ort. 43.4 ± 1.2) |
| hafif_hale | 34.3 | 30.2 | 29.5 | −0.39 /100 Myr (ort. 29.5 ± 0.8) |
| soguk | 44.6 | 35.5 | 37.9 | ~200 Myr periyotla 35–46 arasında salınıyor |

- hafif_hale'deki çubuk sürekli ve düzgün yavaşlamaktadır (34 → 29.5). Bu, kalıcı bir çubuğun dış
  diske açısal momentum verirken beklenen davranıştır.
- Ana koşuda yavaşlama 500–1300 Myr arasında, çubuk kararlıyken gerçekleşmektedir.
- soguk koşusundaki salınım muhtemelen kısa çubukla sarmal kollar arasındaki girişimden
  kaynaklanmaktadır; bu koşunun Ω_p değerleri ve eşdönme oranı güvenilir kabul edilmemiştir.

**Eşdönme / çubuk uzunluğu: ℛ = R_CR / L_çubuk** (t ≥ 500 Myr, çubuk varken):

| Koşu | ℛ, faz-sabit L: medyan [çeyrekler] | ℛ, A2 yarı-maks L | 1.0–1.4 içindeki anlar | n |
|---|---|---|---|---|
| **hafif_hale** | **1.11** [0.98–1.28] | **1.06** [0.94–1.18] | %58 / %55 | 146 |
| ana | 1.41 [1.18–1.66] | 1.60 [1.19–1.95] | %38 / %34 | 88 |
| soguk | 2.10 [1.52–2.45] | 1.42 [1.19–2.15] | %14 / %37 | 146 |

- **hafif_hale gözlenen aralıktadır (1.0–1.4).** İki uzunluk tanımı aynı sonucu vermektedir; çubuk
  "hızlı" sınıfındadır. Güçlü ve kalıcı çubuk, gözlenen çubuklarla uyumlu bir yapıdadır.
- **ana** aralığın üst sınırında ya da biraz dışındadır; zayıf ve kısa ömürlü bir çubukla tutarlıdır.
- **soguk** için iki tanım uyuşmamakta ve Ω_p salınmaktadır; sonuç çıkarılmamıştır.
- Faz-sabit uzunluk, hizalanan sarmal kollar nedeniyle zaman zaman büyüdüğünden medyanlar kullanılmış
  ve ikinci bir tanımla karşılaştırılmıştır.

### 3.4 Merkezi kütle yoğunlaşması

Grafik: `merkezi_kutle.png`.

| Koşu | t = 0 | 500 Myr | 1300 Myr | 1800 Myr | 2000 Myr | 1300→1800 |
|---|---|---|---|---|---|---|
| ana | 2.82 | 4.02 | 5.42 | 5.64 | 5.80 ×10⁹ M☉ | **+%4.0** |
| hafif_hale | 2.82 | 5.16 | 6.58 | 6.91 | 7.05 | +%5.0 |
| soguk | 2.82 | 6.51 | 7.10 | 7.16 | 7.11 | +%0.8 |

Çubuklar merkeze kütle taşımaktadır: merkezi disk kütlesi 2 Gyr'de 2–2.5 katına çıkmaktadır.

---

## 4. Bilinen sınırlar

- **Hale sabit potansiyel olarak modellenmiştir.** Bu nedenle çubuk, dinamik sürtünmeyle haleye
  açısal momentum aktaramaz. Canlı haleli çalışmalarda çubuklar daha güçlü ve uzun ömürlü olmaktadır.
  Ana koşudaki çubuğun zayıflaması da bu sınırdan etkilenmiş olabilir; canlı hale yeni koşu
  gerektirdiğinden test edilmemiştir.
- **Disk ısınması kısmen sayısaldır.** Kalınlık N = 40000'de bile N'ye duyarlıdır (§2.4). Kalınlıkla
  ilgili sayılar, N = 80000 ile bir kontrol yapılmadan fiziksel sonuç olarak yorumlanmamalıdır.
- **hafif_hale'de disk kütlesi yeniden çözülmemiştir.** Bu modelde v_c(8.2) 229 km/s'den düşüktür;
  koşu, yalnızca hale kütlesinin değiştiği bir karşılaştırma olarak tasarlanmıştır.
- **Arayüzdeki analitik Q(R) paneli** Miyamoto–Nagai diskinin M_d'sini (6.8e10) üstel Σ ile
  kullanmaktadır; N-cisim diski 4.82e10'dur.

---

## 5. Açık sorular

### 5.1 Ana koşuda çubuğun zayıflaması

Ana koşudaki çubuğun tarihçesi şöyledir: 480 Myr'de tepe (0.35); 600–1300 Myr arasında ~0.25 civarında
kararlı; **asıl zayıflama 1300–1800 Myr arasında, 0.24 → 0.13.** Bu zayıflama için dört hipotez test
edilmiş ve dördü de elenmiştir. Kontrol olarak çubuğu kalıcı olan hafif_hale koşusu kullanılmıştır.
Grafik: `cubuk_zayiflama_testi.png`.

**a) Bükülme kararsızlığı — elendi.**
- Merkezi bölgenin (R < 3 kpc) kalınlığı tüm koşu boyunca düzgün artmaktadır (162 → 242 pc). Sıçrama
  yoktur: düzgün eğilimden en büyük sapma 11 pc'dir ve t = 10 Myr'dedir.
- Asıl düşüş döneminde (1300–1800 Myr) iç kalınlık %10 artmaktadır (217 → 238 pc); bu, koşunun genel
  eğilimiyle aynıdır.
- A2'nin düştüğü anlarda ortalama kalınlaşma hızı 0.0038 kpc/100 Myr, diğer anlarda 0.0041'dir.
  −dA2/dt ile dz/dt arasındaki korelasyon +0.09'dur.
- İlk küçük düşüşte (480 → 600 Myr) 560–570 Myr'de ~7 pc'lik hafif bir hızlanma vardır (2–3σ). Bu,
  bükülmenin tipik imzası olan, iç kalınlığın ~100 Myr içinde 1.5–2 kat artması değildir.

**b) Disk ısınması — elendi.**

| t [Myr] | Çubuk gücü (50 Myr ort.) | Q (2–4 kpc) | Q (4–8 kpc) |
|---|---|---|---|
| 0 | 0.03 | 1.28 | 1.28 |
| 500 | 0.33 | 1.45 | 1.45 |
| 1200 | 0.25 | 1.74 | 1.76 |
| 1500 | 0.21 | 1.79 | 1.75 |
| 1800 | 0.13 | 1.79 | 1.79 |
| 2000 | 0.16 | 1.81 | 1.82 |

- Q'nun artışının büyük kısmı 0–1200 Myr arasında, **çubuk kararlıyken** gerçekleşmektedir.
- Çubuğun yarıya indiği 1300–1800 Myr arasında Q yalnızca %3 artmaktadır.
- Kontrol koşusu hafif_hale'de Q daha fazla artmakta (2.2'ye), ancak çubuk güçlü kalmaktadır (~0.35).
- Çubuk zayıflarken Q çok az yükselmektedir; ancak ısınma zayıflamanın zamanlamasını açıklamamaktadır.

**c) Desen hızı kaybı — elendi.** 1300–1800 Myr arasında Ω_p düşmemekte, sabit kalmakta ya da çok
hafif artmaktadır (+0.56 /100 Myr; §3.3). Çubuk bu dönemde yavaşlayarak açısal momentum kaybetmemektedir.

**d) Merkezi kütle birikmesi — elendi.**
- Ana koşunun zayıflama döneminde R < 1 kpc içindeki disk kütlesi yalnızca %4 artmaktadır; büyüme hızı
  bu dönemde ortalamanın yarısına düşmektedir (4.6e7'ye karşı 1.1e8 M☉/100 Myr).
- Kontrol koşusu hafif_hale'de merkezde daha fazla kütle birikmekte (7.0e9), ancak çubuk güçlü
  kalmaktadır.
- Ana koşuda M_merkez ile A2 arasındaki −0.84'lük korelasyon, aynı dönemde biri artan biri azalan iki
  eğilimden kaynaklanmaktadır; nedensellik göstermemektedir.
- Ölçülen yoğunlaşma 1 kpc ölçeğindedir. Literatürde çubuğu zayıflatan merkezi kütleler çok daha
  kompakttır (≲100 pc); bu modelde böyle bir yapı yoktur.

**Özet:**

| Hipotez | 1300–1800 Myr arasında gözlenen | Sonuç |
|---|---|---|
| Bükülme | İç kalınlık düzgün artıyor, sıçrama yok | elendi |
| Disk ısınması | Q yalnızca %3 artıyor; asıl artış çubuk kararlıyken oluyor | elendi |
| Desen hızı kaybı | Ω_p düşmüyor (+0.56 /100 Myr) | elendi |
| Merkezi kütle | M(<1 kpc) yalnızca %4 artıyor; kontrolde daha fazla kütle var, çubuk kalıcı | elendi |

**Zayıflamanın nedeni açık bir sorudur.** Test edilmemiş olasılıklar: sabit hale (§4), sayısal ısınma
(N'ye bağlılık) ve çubuk–sarmal etkileşimi. Üçü de yeni koşu gerektirmektedir: canlı haleli bir koşu,
N = 80000'lik bir koşu ve sarmal yapının bastırıldığı bir kontrol.

### 5.2 Disk kalınlığının yakınsaması

Kalınlık N = 40000'de tam yakınsamamıştır (§2.4). N = 80000 ile bir kontrol, kalınlaşmanın fiziksel
ve sayısal kısımlarını ayırmak için gereklidir.

---

## 6. Grafikler ve veriler

Depoda bulunanlar (`runs/`):

- `a2_karsilastirma.png` — çubuk gücü ve uzunluğu
- `cubuk_faz_profilleri.png` — A2/A0(R) ve m=2 fazı(R)
- `cubuk_zayiflama_testi.png` — bükülme ve ısınma testleri
- `ana/animasyon.mp4`, `hafif_hale/animasyon.mp4`

Koşulardan yeniden üretilebilenler (`.venv/bin/python -m src.analysis`, ~40 s):

- `desen_hizi.png` — Ω_p(t) ve R_CR / L_çubuk
- `merkezi_kutle.png` — M(R < 1 kpc)(t) ile çubuk gücü
- `a2_toplu_R5.png`, `dogrulama_karsilastirma.png`
- `<ad>/karsilastirma.png` — 0/500/1000/2000 Myr; doğrulama koşularında 0/200/400/600 Myr
- `soguk/Q_profili.png`, `soguk/parcalanma.png`
- `ozet_metrikler.json` — tüm sayısal metrikler

Küçültülmüş koşu verisi (3B görüntüleme için, 5000 parçacık, 80 Myr aralık): `deploy/runs/`.
