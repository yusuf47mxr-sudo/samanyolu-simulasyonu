# Samanyolu Simülasyonu — Proje Tanımı

Bu belge projenin fizik modelini, parametrelerini, birim sistemini, kabul kriterlerini ve
sayısal kararlarını tanımlar. Sonuçlar ve ayrıntılı doğrulama [runs/SONUCLAR.md](runs/SONUCLAR.md)
belgesindedir.

---

## 1. Amaç

Proje, Samanyolu'nun dinamiğini iki aşamada simüle eder ve fiziksel parametreleri canlı
değiştirilebilen bir arayüzle sunar:

1. **Test parçacık modeli** — sabit analitik potansiyelde yıldız yörüngeleri ve dönme eğrisi.
2. **N-cisim modeli** — öz-yerçekimli bir disk; çubuk ve sarmal kolların kendiliğinden oluşumu.

Projenin bilimsel mesajı: **galaksiyi bir arada tutan şey görünür madde değildir.** Arayüz bu
mesajı göstermek üzere tasarlanmıştır.

---

## 2. Birim sistemi

Tüm iç hesaplar aşağıdaki birimlerde yapılır; dönüşüm yalnızca ekrana yazılırken uygulanır.

| Büyüklük | Birim |
|---|---|
| Uzunluk | kpc |
| Zaman | Myr |
| Kütle | M☉ |

```python
G        = 4.498502e-12    # kpc^3 / (M_sun * Myr^2)
KMS_KPMY = 1.022712e-3     # 1 km/s  ->  kpc/Myr
```

Hızlar içeride kpc/Myr cinsindendir; km/s'ye yalnızca grafik etiketlerinde ve arayüzde çevrilir.

---

## 3. Fizik modeli

Toplam potansiyel dört bileşenin toplamıdır; kuvvetler de bileşen bazında toplanır.

### 3.1 Disk — Miyamoto–Nagai

Bu projede disk Miyamoto–Nagai potansiyeli ile modellenmiştir:

```
Phi(R,z) = -G*M_d / sqrt( R^2 + (a + sqrt(z^2 + b^2))^2 )
```

| Parametre | Varsayılan | Aralık |
|---|---|---|
| `M_d` | 6.8e10 M☉ | 3e10 – 12e10 |
| `a` (ölçek uzunluğu) | 3.0 kpc | 1.5 – 6.0 |
| `b` (ölçek yüksekliği) | 0.28 kpc | 0.1 – 1.0 |

### 3.2 Şişkinlik — Hernquist

```
Phi(r) = -G*M_b / (r + c)
```

| Parametre | Varsayılan | Aralık |
|---|---|---|
| `M_b` | 5.0e9 M☉ | 0 – 2e10 |
| `c` | 0.5 kpc | 0.1 – 2.0 |

### 3.3 Karanlık madde halesi — NFW

```
rho(r)   = rho_s / [ (r/r_s) * (1 + r/r_s)^2 ]
M(<r)    = 4*pi*rho_s*r_s^3 * [ ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s) ]
Phi(r)   = -4*pi*G*rho_s*r_s^3 * ln(1 + r/r_s) / r
```

Normalizasyon (`conc = r_vir / r_s`):

```
rho_s * r_s^3 = M_vir / ( 4*pi * [ ln(1+conc) - conc/(1+conc) ] )
```

| Parametre | Varsayılan | Aralık |
|---|---|---|
| `M_vir` | 1.0e12 M☉ | 0 – 3e12 |
| `r_s` | 16.0 kpc | 8 – 30 |
| `conc` | 15.0 | 5 – 25 |

Hale bir aç/kapa anahtarına sahiptir; kapalıyken bu bileşen potansiyele hiç katılmaz.

### 3.4 Merkezi kara delik — Kepler

```
Phi(r) = -G*M_bh / r
```

| Parametre | Varsayılan | Aralık |
|---|---|---|
| `M_bh` | 4.3e6 M☉ | 0 – 1e10 |

Üst sınır bilerek abartılmıştır. Gerçek değerde kara deliğin dönme eğrisine etkisi görünmez; bu,
modelin gösterdiği bir sonuçtur. Etki yarıçapı (~2 pc) yumuşatma uzunluğundan küçüktür. N-cisim
koşularında ve canlı modda kara delik, parçacıklarla aynı uzunlukla (eps) Plummer yumuşatmasına tabidir.

### 3.5 Dönme eğrisi

```
v_c^2(R) = R * dPhi/dR  |_(z=0)      (bileşenler toplanarak)
```

Türevler analitiktir. Sayısal türev yalnızca testlerde karşılaştırma amacıyla kullanılmıştır.

---

## 4. Kabul kriterleri

Her kriter `tests/` altında otomatik bir testtir.

| # | Test | Beklenen |
|---|---|---|
| T1 | Varsayılan parametrelerle `v_c(8.2 kpc)` | **229 ± 5 km/s** |
| T2 | Dönme eğrisi 5–20 kpc arası | 205 – 240 km/s bandında, düz |
| T3 | Güneş'in dairesel yörünge periyodu `2*pi*R/v_c` | **220 ± 10 Myr** |
| T4 | Eksenel simetrik potansiyelde tek yörünge, `L_z` korunumu | bağıl sapma < 1e-6 |
| T5 | Leapfrog ile 1000 adımda toplam enerji | bağıl sapma < 1e-3 |
| T6 | Hale kapalıyken `v_c(20 kpc)` | **< 165 km/s** (Kepler düşüşü görünür) |
| T7 | `M_bh`'yi 4.3e6 → 4.3e8 yap, `v_c(8.2)` değişimi | < 1 km/s |
| T8a | Dinamik kütle `M_dyn = R*v_c(R)^2/G`, R = 8.2 kpc | **1.00e11 ± 0.05e11 M☉** (gözlem makaleleriyle karşılaştırılabilir olan) |
| T8b | Küre içi gerçek kütle: `rho`'nun r < 8.2 kpc küresi üzerinde integrali | **0.90e11 ± 0.05e11 M☉** |
| T9 | N-cisim diskinin kendi parçacıklarından ölçülen dönme hızı, 8.2 kpc | **229 ± 5 km/s** |

**T8 notu.** `M_dyn / M_true ≈ 1.11` çıkar; fark diskin yassılığından gelir. Küresel bir kütle
dağılımında kabuk teoremi gereği oran tam 1.00 olurdu. T8a gözlemle karşılaştırılacak olandır.
T8b, `rho`'nun küre üzerindeki integralidir; disk için Gauss yasasıyla hesaplanır ve doğrudan hacim
integraliyle aynı sonucu verdiği ayrıca test edilir.

**T9 notu.** v_c, parçacıkların ve sabit bileşenlerin 8.2 kpc halkasındaki kütleçekiminden ölçülür;
ortalama ⟨v_phi⟩ değildir. ⟨v_phi⟩ asimetrik sürüklenme kadar düşüktür (~9 km/s) ve bu fark Jeans
denkleminin öngörüsüyle %20 içinde uyuşur (testli).

T5'in uzun koşu (2000 Myr, N-cisim) hâli de bir kabul kriteri olarak izlenmiştir.

---

## 5. Proje aşamaları

Proje dört aşamada geliştirilmiştir. Her aşama kendi kabul testleri geçtikten sonra kapatılmıştır.

### Aşama 0 — Potansiyeller
`src/units.py`, `src/potentials.py`, `src/rotation.py`. Analitik potansiyeller, dönme eğrisi ve
bileşen bazında dönme eğrisi grafiği. Testler: T1, T2, T6, T7, T8a, T8b.

### Aşama 1 — Yörüngeler ve arayüz
`src/orbits.py` (leapfrog), `app.py` (Streamlit). Çalışan, kaydıraklı arayüz. Testler: T3, T4, T5.
Bu aşama tek başına kullanılabilir bir üründür.

### Aşama 2 — N-cisim
`src/ics.py`, `src/nbody.py`. Hale **sabit potansiyel** olarak kalmış, canlı parçacık
kullanılmamıştır; disk öz-yerçekimlidir. Testler: T5 (uzun koşuda), T9 ve Q < 1.0 ile başlatılan
diskin parçalanması.

### Aşama 3 — Analiz
`src/analysis.py`. Çubuk gücü A2/A0'ın zamanla değişimi, desen hızı ve eşdönme, animasyon üretimi
ve önceden koşulmuş senaryoların arayüze bağlanması.

**Çubuk ölçümü.** A2/A0, 0.5 kpc'lik kuşaklarda hesaplanır; çubuk gücü bunun R < 5 kpc'deki
maksimumudur. Çubuk uzunluğu, m=2 fazının ±10° sabit kaldığı bölgedir. Bir yapının çubuk sayılması
için güç ≥ 0.2 ve faz-sabit uzunluk ≥ 2 kpc gerekir; böylece sarılmış sarmallar çubuk sayılmaz.
Oluşum zamanı, bu ölçütün 100 Myr kesintisiz sağlandığı ilk andır.

---

## 6. Başlangıç koşulları (N-cisim)

- Yüzey yoğunluğu: `Sigma(R) = Sigma_0 * exp(-R/R_d)`, `R_d = 2.6 kpc`, 15 kpc'de kesik.
- Dikey dağılım: `sech^2(z/z_0)`, `z_0 = 0.3 kpc`.
- Dairesel hız toplam potansiyelden alınır; asimetrik sürüklenme düzeltmesi uygulanır.
- Radyal hız dağılımı hedef Toomre Q'dan türetilir:

```
Q = sigma_R * kappa / (3.36 * G * Sigma)
kappa^2 = R * d(Omega^2)/dR + 4*Omega^2
```

| Parametre | Varsayılan | Aralık | Not |
|---|---|---|---|
| `Q_target` | 1.3 | 0.8 – 2.5 | <1 parçalanır, >2 kol oluşmaz |
| `N` | 20000 | 5000 – 50000 | sayısal |
| `eps` (yumuşatma) | 0.08 kpc | 0.05 – 0.15 | sayısal |
| `dt` | 0.2 Myr | 0.05 – 0.5 | sayısal |
| `t_end` | 2000 Myr | 500 – 3000 | fiziksel |
| `theta` (Barnes–Hut) | 0.6 | 0.4 – 0.8 | sayısal |

Disk/hale kütle oranı çubuk oluşumunu belirleyen ana parametredir (Ostriker–Peebles).

**Disk kütlesinin yeniden normalizasyonu.** `R_d = 2.6 kpc` ölçülmüş bir büyüklüktür; `M_d = 6.8e10`
ise Miyamoto–Nagai uyum fonksiyonunun parametresiydi. Fonksiyon biçimi üstele değiştiği için
normalizasyon yeniden yapılmıştır; korunan büyüklük gözlemdir: `v_c(8.2) = 229 km/s`. Üstel disk
kütlesi sabit yazılmamış, `ics.solve_disk_mass` ile (brentq) üstel disk + şişkinlik + hale toplamı
8.2 kpc'de 228.7 km/s verecek şekilde çözülmüştür: **M_d ≈ 4.82e10 M☉**. Bu çözüm yalnızca varsayılanı
belirler (`config.DiskICParams.M_d`).

**Çubuk kararlılığı tanısı.** Her başlangıç koşulu üretiminde Efstathiou–Lake–Negroponte parametresi
hesaplanır: `eps = V_max / sqrt(G*M_d/R_d)`. `eps < 1.1` çubuk oluşumuna işaret eder; varsayılan
model için eps ≈ 0.81.

**Başlangıçta dairesel hız.** Disk kalın (`sech^2`) ve yumuşatılmış olduğundan düzlemdeki v_c, ince
(Freeman) diskinkinden düşüktür (iç diskte 8 km/s'ye kadar). Bu nedenle ⟨v_phi⟩, parçacıkların
kendi alanından halkalarda ölçülen v_c ile kurulmuştur.

### 6.1 Sayısal kararlar

Aşağıdaki kararlar bir yakınsama çalışmasıyla belirlenmiştir (ayrıntılar SONUCLAR.md'de).

- **Barnes–Hut ağacı kuadrupol momentlidir**, `theta = 0.6`. Monopol ağaç enerjiyi tek yönde
  sürüklüyordu (600 Myr'de 4e-4; 2000 Myr'de T5 aşılırdı). Kuadrupol ile RMS kuvvet hatası ~%0.1'e
  (hedef < %0.3), enerji sapması 2000 Myr'de 3e-5 – 8e-5'e inmiştir.
- **Üretim koşuları N = 40000** ile yapılmıştır; varsayılan 20000 hızlı denemeler içindir. θ ve dt
  sonucu değiştirmemiş, N disk kalınlığını değiştirmiştir (iki-cisim ısınması). N = 40000'de bile
  kalınlık tam yakınsamış değildir; kalınlıkla ilgili sonuçlar N = 80000 ile kontrol edilmeden
  fiziksel kabul edilmez.
- `eps = 0.08 kpc` ve `dt = 0.2 Myr` değişmemiştir.
- Yaklaşık yöntem (ağaç) tam yöntemle (doğrudan N²) karşılaştırılmıştır: R_d ve A2/A0 farkları
  yalnızca seed değişince oluşan farktan küçüktür.

### 6.2 Modelin bilinen sınırı: sabit hale

Hale sabit potansiyeldir. Bu nedenle çubuk, dinamik sürtünmeyle haleye açısal momentum aktaramaz.
Canlı haleli çalışmalarda çubuklar daha güçlü ve uzun ömürlü olmaktadır. Sonuçlar bu sınırla
birlikte yorumlanmalıdır.

---

## 7. Arayüz

Streamlit arayüzünde sol sütunda kaydıraklar, sağda canlı güncellenen paneller bulunur.

**Kontroller.** Ana sekmede beş kontrol vardır: `M_d`, `M_vir`, karanlık madde aç/kapa, `Q_target`,
`M_bh`. "Gelişmiş" sekmesi diğer fiziksel parametreleri, "Sayısal" sekmesi `N`, `eps`, `dt`, `theta`
parametrelerini içerir. Sayısal sekmenin başında şu not yer alır: *bu parametreler sonucu
değiştirmemelidir; değiştiriyorsa simülasyon yakınsamamıştır.*

**Analitik model sekmesi:**
1. Dönme eğrisi — toplam ve bileşenler ayrı çizgilerle, üzerinde gözlem noktaları.
2. Logaritmik yarıçap paneli (0.001 – 300 kpc) — kara delik, şişkinlik, disk ve hale rejimlerinin
   dördü birden görünür.
3. Sayısal özet: `v_c(8.2)`, yörünge periyodu, `M(<8.2 kpc)`, kaçış hızı.
4. `Q(R)` profili.

**Ters mod.** Kullanıcı hedef `v_c(8.2)` girer, gerekli `M_vir` `scipy.optimize.brentq` ile çözülür.

**Canlı sekmesi.** Test parçacıkları seçili potansiyelde gerçek zamanlı ve üç boyutlu olarak
ilerletilir; kaydıraklar anında etkir. Öz-yerçekimi yoktur.

**Koşu sonuçları sekmesi.** Önceden koşulmuş N-cisim simülasyonlarının 3B görünümü, zaman kaydırağı,
A2/A0 zaman serisi ve enerji korunum grafiği.

---

## 8. Tasarım ilkeleri

- **Zaman entegrasyonu leapfrog (kick–drift–kick) ile yapılır.** RK4 simplektik olmadığından enerji
  sürükler; kullanılmamıştır.
- **Koşu sırasında çizim yapılmaz.** Belirli adım aralıklarında `runs/<ad>/snap_XXXX.npy` yazılır;
  animasyon koşu bittikten sonra bu dosyalardan üretilir.
- **Fiziksel ve sayısal parametreler ayrı tutulur:** ayrı dataclass, ayrı arayüz sekmesi, ayrı belge
  bölümü.
- Kuvvet hesabı döngüsü `numba.njit` ile derlenir.
- Tüm varsayılanlar `src/config.py` içinde tek yerde durur.
- Her fonksiyonun belge dizisinde birimler yazılıdır (`R : kpc`, `v : kpc/Myr`).
- Simülasyon motoru düz `.py` modüllerinden oluşur.

---

## 9. Dosya yapısı

```
MilkyWay/
├── PROJE_TANIMI.md         # bu belge
├── README.md
├── LICENSE
├── requirements.txt
├── app.py                  # Streamlit arayüzü
├── src/
│   ├── units.py            # sabitler, dönüşümler
│   ├── config.py           # tüm varsayılan parametreler
│   ├── potentials.py       # MN, Hernquist, NFW, Kepler
│   ├── rotation.py         # dönme eğrisi, kappa, Omega
│   ├── orbits.py           # leapfrog, tek parçacık yörüngeleri
│   ├── ics.py              # N-cisim başlangıç koşulları
│   ├── nbody.py            # kuvvet hesabı, zaman adımı
│   └── analysis.py         # A2/A0, enerji, animasyon
├── ui/galaxy3d/            # 3B görünüm bileşeni (kamera açısını korur)
├── deploy/                 # bulut dağıtımı: app.py giriş noktası, runs/ küçültülmüş veri
├── docs/                   # ekran görüntüleri
├── tests/
│   ├── test_potentials.py  # T1, T2, T6, T7, T8a, T8b
│   ├── test_orbits.py      # T3, T4
│   ├── test_energy.py      # T5
│   ├── test_ics.py         # T9, Q, asimetrik sürüklenme
│   ├── test_nbody.py       # N-cisim enerji (T5)
│   ├── test_tree_accuracy.py  # ağaç kuvvet doğruluğu
│   ├── test_analysis.py    # çubuk ölçümleri
│   └── test_deploy.py      # küçültülmüş veri: boyut < 20 MB, biçim
└── runs/                   # snapshot çıktıları (depoda yalnızca SONUCLAR.md ve seçili grafikler)
```

Bağımlılıklar (`requirements.txt`, sürümleri sabit): `numpy scipy matplotlib numba streamlit plotly pytest`.

---

## 10. Geliştirme yöntemi

Proje aşamalı geliştirilmiştir: her aşama testleri çalıştırılıp sonuçları incelendikten sonra bir
sonrakine geçilmiştir. Bir kabul kriteri tutmadığında kriter gevşetilmemiş, kod düzeltilmiştir.
Kriterin kendisinin yanlış ya da eksik tanımlandığı durumlarda (ör. T8'in T8a/T8b olarak ayrılması)
değişiklik gerekçesiyle birlikte bu belgeye işlenmiştir.
