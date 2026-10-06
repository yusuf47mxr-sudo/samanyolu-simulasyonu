# Aşama 2c + Aşama 3 — Özet

Son güncelleme: 2026-10-04 (desen hızı, eşdönme ve merkezi kütle eklendi). Bütün testler:
**50/50 geçti** (`.venv/bin/python -m pytest tests`). Yeni koşu yapılmadı; bütün sonuçlar mevcut
snapshot'lardan. **Fizik tarafı bu güncellemeyle kapandı.**

**Öne çıkan sonuç:** hafif_hale çubuğunun eşdönme / uzunluk oranı **ℛ = 1.11** (faz-sabit uzunlukla;
yarı-maksimum uzunlukla 1.06). Bu, gözlenen çubukların **1.0–1.4** aralığının içinde: hızlı bir çubuk
(§2.6).

## Sabah bakılacaklar (güncel)

1. **Çubuk ölçümü yeniden tanımlandı** (karar 1): A2/A0 artık radyal kuşaklarda ölçülüyor ve
   çubuk uzunluğu m=2 fazının sabit kaldığı bölgeden geliyor. Bunu uygularken tanıma iki ek koydum;
   **onayınız gerekiyor** (bkz. §2.1):
   - Maksimum, toplu ölçümle aynı R < 5 kpc bölgesindeki kuşaklarda alınıyor. Tüm yarıçaplardaki
     maksimum sarmal kolları yakalıyor; o da grafikte ince çizgiyle gösteriliyor.
   - Bir m=2 yapısının "çubuk" sayılması için fazın en az 2 kpc boyunca sabit kalması gerekiyor.
     Bu koşul olmadan merkezdeki sarılmış sarmallar 50–120 Myr'de sahte "çubuk oluşumu" veriyordu.
2. **Soğuk disk: kriter karşılandı, öz-düzenleme ölçüldü** (§2.2): Q, 0.80'den 800 Myr'de 2.32'ye çıkıyor.
3. **Ana koşuda çubuğun zayıflaması: açıklanamadı** (§2.3, §2.6, §2.7). Beş aday test edildi:
   bükülme, ısınma, desen hızı kaybı, merkezi kütle birikmesi ve (dolaylı olarak) çubuk uzunluğu ölçümü.
   Hiçbiri asıl düşüş dönemini (1300–1800 Myr) açıklamıyor.
4. **Model sınırı** (§3): hale sabit potansiyel.

---

## 1. Ne yapıldı

| Adım | Durum |
|---|---|
| RUN-1 `ana` (N=40000, θ=0.6+kuadrupol, Q=1.3, 2000 Myr) | ✓ 1488 s |
| RUN-2 `soguk` (Q=0.8) | ✓ 1613 s |
| RUN-3 `hafif_hale` (M_vir = 5e11) | ✓ 1508 s |
| RUN-4 `dogrulama_dogrudan` / `dogrulama_agac` (N=5000, 600 Myr, aynı seed) | ✓ |
| + `dogrulama_dogrudan_seed7` (gürültü tabanı için, benim eklediğim) | ✓ |
| `src/analysis.py`: A2/A0(R), m=2 fazı, çubuk uzunluğu, Q(R), enerji, Σ(R), kalınlık, animasyon | ✓ |
| `tests/test_analysis.py` (sentetik veri + koşu tutarlılığı, 15 test) | ✓ |
| Grafikler ve animasyonlar (yeni çubuk ölçüsüyle yeniden üretildi) | ✓ |
| `app.py` → "Koşu sonuçları" sekmesi (A2/A0(R) ve faz profili dahil) | ✓ |

**Dokunulmayanlar:** `potentials.py`, `ics.py`, `nbody.py`, `config.py`, `CLAUDE.md`. Hiçbir koşu
klasörü silinmedi. Her koşu klasörüne yalnızca türetilmiş dosyalar eklendi: `bar_history.npz` ve `pattern_history.npz`
(önbellek), `karsilastirma.png`, ana ve hafif_hale'ye `animasyon.mp4`, soğuğa `parcalanma.png` ve
`Q_profili.png`.

**Koşuların bittiği nasıl anlaşıldı:** Doğru sayı 2000/10 + 1 = **201 snapshot** (+ 4 yardımcı
dosya). `analysis.is_complete()` bunu config.json'dan hesaplıyor.

## 2. Sonuçlar

### 2.1 Çubuk: güç, uzunluk, oluşum zamanı

Grafikler: `runs/a2_karsilastirma.png` (üstte güç, altta uzunluk), `runs/cubuk_faz_profilleri.png`
(A2/A0(R) ve faz(R)).

**Yöntem:**
- A2/A0, 0.5 kpc'lik kuşaklarda hesaplanıyor.
- **Çubuk gücü** = A2/A0(R)'nin R < 5 kpc'deki maksimumu. Yalnızca en az 200 parçacıklı kuşaklar
  sayılıyor.
- **Çubuk uzunluğu** = tepe kuşağından dışarı doğru m=2 fazının tepe fazından ±10° içinde kaldığı son
  kuşağın dış kenarı.
- **Çubuk var** = güç ≥ 0.2 ve uzunluk ≥ 2 kpc.
- **Oluşum zamanı** = bu koşulun 100 Myr kesintisiz sağlandığı ilk an.
- Gürültü tabanı: t = 0'da güç 0.034 (N=40k), 0.06–0.12 (N=5k).

| Koşu | Oluşum | Çubuk var (tüm / son %25) | Tepe gücü | Son %25 ort. güç | Uzunluk medyanı (son %25) |
|---|---|---|---|---|---|
| ana | 390 Myr | %54 / **%14** | 0.38 @ 480 Myr | 0.15 | 3.5 (3.0) kpc |
| soguk | **140 Myr** | %98 / %100 | 0.50 @ 250 Myr | 0.24 | 2.5 (2.5) kpc |
| hafif_hale | **150 Myr** | %91 / **%100** | **0.50** @ 480 Myr | **0.35** | **5.0 (6.0) kpc** |

**Beklenti tuttu:** hafif_hale'de çubuk ana koşudan **240 Myr daha erken** oluşuyor, **daha güçlü**
ve **daha uzun**: uzunluk medyanı 5.0 kpc, ana koşuda 3.5 kpc (~1.4 kat). 2000 Myr boyunca da kalıcı. Ana koşuda çubuk oluşuyor ama son
çeyrekte büyük ölçüde sönmüş.

Notlar:
- **Uzunluk ölçüsü anlık olarak gürültülü** (2–9 kpc arası sıçrıyor): fazı çubukla hizalanan bir
  sarmal kol, çubuğu "uzatıyor". Medyanlar güvenilir, tekil değerler değil.
- **N=5000 doğrulama koşularında oluşum zamanı kaotik:** aynı seed'le doğrudan toplamada 160 Myr,
  ağaçta 390 Myr; farklı seed'le doğrudan toplamada 180 Myr. Bu, N=5000'de oluşum zamanının
  yöntemden çok gürültüye bağlı olduğunu gösteriyor. Ortalama güç ve uzunluk için §2.4'e bakın.
- Eski tanım (R < 5 kpc toplu A2) karşılaştırma için `runs/a2_toplu_R5.png`'de duruyor.

### 2.2 Soğuk disk (Q = 0.8): parçalanma ve öz-düzenleme

Grafikler: `runs/soguk/Q_profili.png`, `runs/soguk/parcalanma.png`.

**Kriter karşılandı.** Disk ilk 200 Myr'de parçalanıyor: topak kesri 150–200 Myr'de 0.20'ye
çıkıyor; ana koşuda aynı anda 0.05–0.06.

**Öz-düzenleme:** Kararsız disk kendi kendini ısıtıyor ve Q kararlı bölgeye yükseliyor.

| t | Ortalama Q (1–12 kpc) | Aralık |
|---|---|---|
| 0 Myr | **0.80** | 0.76–0.84 |
| 200 Myr | 1.80 | 1.50–2.24 |
| 800 Myr | **2.32** | 2.05–2.92 |

800 Myr'de **her yarıçapta Q ≥ 2.0**. Topaklar 500 Myr'de kaybolmuş, disk yayılmış (R_d 2.6 → 4.3 kpc).
Çarpışmasız yıldız diski kararsızlığı ısınarak doyuruyor; gaz diskinden farklı olarak kalıcı topak
oluşmuyor.

**Q nasıl ölçüldü:** Σ ve σ_R parçacıklardan, κ ise parçacıkların ve sabit bileşenlerin halkalar
üzerinde ölçülen v_c'sinden hesaplandı. Disk zamanla yayıldığı için başlangıçtaki analitik κ
kullanılmadı. t=0'da bu yöntem analitik κ ile R > 4 kpc'de %1 içinde uyuşuyor; içeride %4–6 düşük
(kalın diskin gerçek alanı). Çubuk bölgesinde σ_R akış hareketlerini de içeriyor; orada Q kabaca
anlamlı.

### 2.3 Ana koşuda çubuğun zayıflaması: iki hipotez

Grafik: `runs/cubuk_zayiflama_testi.png`. Satırlar ana ve hafif_hale (kontrol); sol sütun bükülme,
sağ sütun ısınma testi.

Yeni ölçüye göre ana koşudaki çubuğun tarihçesi:
- 480 Myr'de tepe (0.35);
- 600–1300 Myr arasında ~0.25 civarında kararlı;
- **asıl zayıflama 1300–1800 Myr arasında: 0.24 → 0.13.**

**a) Bükülme kararsızlığı: TUTMUYOR.**
- Merkezi bölgenin (R < 3 kpc) kalınlığı tüm koşu boyunca düzgün artıyor (162 → 242 pc). Sıçrama yok:
  düzgün eğilimden en büyük sapma 11 pc ve o da t = 10 Myr'de.
- Asıl düşüş döneminde (1300–1800 Myr) iç kalınlık %10 artıyor (217 → 238 pc); bu, koşunun genel
  eğilimiyle aynı.
- A2'nin düştüğü anlarda ortalama kalınlaşma hızı 0.0038 kpc/100 Myr, diğer anlarda 0.0041.
  −dA2/dt ile dz/dt arasındaki korelasyon +0.09.
- İlk küçük düşüşte (480 → 600 Myr) 560–570 Myr'de ~7 pc'lik hafif bir hızlanma var (2–3σ). Bükülmenin
  tipik imzası olan, iç kalınlığın ~100 Myr içinde 1.5–2 kat artması değil.

**b) Disk ısınması: asıl düşüşü AÇIKLAMIYOR.**

| t [Myr] | Çubuk gücü (50 Myr ort.) | Q (2–4 kpc) | Q (4–8 kpc) |
|---|---|---|---|
| 0 | 0.03 | 1.28 | 1.28 |
| 500 | 0.33 | 1.45 | 1.45 |
| 1200 | 0.25 | 1.74 | 1.76 |
| 1500 | 0.21 | 1.79 | 1.75 |
| 1800 | 0.13 | 1.79 | 1.79 |
| 2000 | 0.16 | 1.81 | 1.82 |

- Q'nun artışının büyük kısmı 0–1200 Myr arasında, **çubuk kararlıyken** oluyor.
- Çubuğun yarıya indiği 1300–1800 Myr arasında Q yalnızca %3 artıyor.
- Kontrol koşusu hafif_hale'de Q daha da fazla artıyor (2.2'ye), ama çubuk güçlü kalıyor (~0.35).
- "Çubuk zayıflarken Q yükseliyor mu?" sorusunun cevabı kelimenin tam anlamıyla evet (çok az), ama ısınma
  zayıflamanın zamanlamasını açıklamıyor.

**Bu iki hipotezin sonucu:** ikisi de asıl zayıflamayı açıklamıyor. Desen hızı ve merkezi kütle
testleri için §2.6 ve §2.7'ye bakın; nihai sonuç §2.7'nin sonunda.

### 2.4 Çapraz doğrulama: ağaç vs doğrudan toplama

Grafik: `runs/dogrulama_karsilastirma.png`.

| Karşılaştırma | Ortalama \|ΔR_d\| | Ortalama \|Δ(A2/A0)\| (toplu) |
|---|---|---|
| ağaç+kuadrupol − doğrudan (aynı seed) | **0.065 kpc** | **0.037** |
| doğrudan seed 42 − doğrudan seed 7 (gürültü tabanı) | 0.211 kpc | 0.041 |

Yaklaşık yöntemle tam yöntem arasındaki fark gürültü tabanından küçük: **"yaklaşık yöntemimiz tam
yöntemle istatistiksel olarak aynı sonucu veriyor"**. Yeni ölçüyle son %25'teki ortalama çubuk gücü
0.25 (ağaç), 0.37 (doğrudan), 0.29 (seed 7); uzunluk medyanları 2.5–3.0 kpc ve üçünde aynı. Güçteki
fark seed farkıyla aynı mertebede.

### 2.5 Enerji ve kütle merkezi

| Koşu | Maks \|ΔE/E\| (2000 Myr) | Maks \|r_com\| | Maks \|v_com\| |
|---|---|---|---|
| ana | **3.0e-5** | 0.15 kpc | 4.1 km/s |
| soguk | **8.1e-5** | 0.27 kpc | 12.7 km/s |
| hafif_hale | **5.8e-5** | 0.24 kpc | 5.3 km/s |

- T5 sınırı (1e-3) uzun koşularda 12–33 kat altında kaldı.
- **Kütle merkezi sürekli kaymıyor, salınıyor.** Soğuktaki 0.27 kpc'lik sapma parçalanma evresinde
  (200–300 Myr) oluşuyor; 500 Myr'den sonra ~0.05 kpc'ye geri dönüyor.

### 2.6 Desen hızı ve eşdönme

Grafik: `runs/desen_hizi.png`.

**Yöntem:**
- **Çubuk fazı:** m=2 fazı, faz-sabit çubuk bölgesinden (0.5 kpc < R < L_çubuk) ölçülüyor.
- **Desen hızı:** Ω_p, açılmış fazın ±50 Myr penceresindeki doğrusal eğimi. Bir m=2 fazı 180°'de bir
  tekrarladığı için açma, adım başına < 90° dönmeyi varsayıyor. t ≥ 500 Myr'de en büyük adım
  23–41°, yani açma güvenli. Daha erken dönemde 60–89°'lik adımlar var; o yüzden Ω_p yalnızca
  t ≥ 500 Myr için veriliyor.
- **Eşdönme yarıçapı:** R_CR, o snapshot'ın gerçek alanından ölçülen v_c ile çözülüyor:
  v_c(R_CR) = Ω_p · R_CR.

**Desen hızı:**

| Koşu | Ω_p, 500 Myr | Ω_p, 1300 Myr | Ω_p, 1800 Myr | 1300–1800 eğilimi |
|---|---|---|---|---|
| ana | 48.6 | 41.4 | 44.1 km/s/kpc | **+0.56 /100 Myr** (ort. 43.4 ± 1.2) |
| hafif_hale | 34.3 | 30.2 | 29.5 | −0.39 /100 Myr (ort. 29.5 ± 0.8) |
| soguk | 44.6 | 35.5 | 37.9 | ~200 Myr periyotla 35–46 arasında salınıyor |

- **ana, 1300–1800 Myr: Ω_p düşmüyor**, sabit ya da çok hafif artıyor. Yani çubuk bu dönemde
  yavaşlayarak açısal momentum kaybetmiyor. Yavaşlama 500–1300 Myr arasında, çubuk kararlıyken oluyor.
- hafif_hale'deki çubuk sürekli ve düzgün yavaşlıyor (34 → 29.5). Bu, kalıcı bir çubuğun dış diske
  açısal momentum verirken beklenen davranış.
- **soguk'taki salınım:** muhtemelen kısa çubukla sarmal kollar arasındaki girişim. Bu koşunun Ω_p
  değerlerini ve eşdönme oranını güvenilir saymıyorum.

**Eşdönme / çubuk uzunluğu: ℛ = R_CR / L_çubuk** (t ≥ 500 Myr, çubuk varken):

| Koşu | ℛ, faz-sabit L: medyan [çeyrekler] | ℛ, A2 yarı-maks L | 1.0–1.4 içindeki anlar | n |
|---|---|---|---|---|
| **hafif_hale** | **1.11** [0.98–1.28] | **1.06** [0.94–1.18] | %58 / %55 | 146 |
| ana | 1.41 [1.18–1.66] | 1.60 [1.19–1.95] | %38 / %34 | 88 |
| soguk | 2.10 [1.52–2.45] | 1.42 [1.19–2.15] | %14 / %37 | 146 |

- **hafif_hale: gözlenen aralıkta (1.0–1.4).** İki uzunluk tanımı da aynı sonucu veriyor; çubuk
  "hızlı" sınıfında. Güçlü ve kalıcı çubuk, gözlenen çubuklarla uyumlu bir yapıda.
- **ana:** aralığın üst sınırında ya da biraz dışında; zayıf ve kısa ömürlü bir çubuk için tutarlı.
- **soguk:** iki tanım uyuşmuyor ve Ω_p salınıyor; sonuç çıkarmıyorum.
- **Uzunluk ölçüsü uyarısı:** Faz-sabit uzunluk, hizalanan sarmal kollar yüzünden zaman zaman
  büyüyor. Bu yüzden medyanlar kullanıldı ve ikinci bir tanımla karşılaştırıldı.

**hafif_hale çubuk uzunluğu kesilmiş mi? Hayır, tablo değişmiyor.**
- 5 kpc sınırı yalnızca tepe aramasında vardı; faz yürüyüşü zaten 10 kpc'ye kadar gidiyordu.
  Çubuk olan anların %46'sında L > 5 kpc ölçülmüş.
- Arama 10 kpc'ye açılınca tepe yalnızca 12 snapshot'ta 5 kpc'nin dışına kayıyor; hepsi t ≤ 310 Myr'deki
  sarmal evresinde.
- Faz çoğu anda 5–5.5 kpc'den sonra düzgünce kayıyor (çubuktan sarmala geçiş). 1700–1800 Myr'de faz
  7 kpc'ye kadar sabit, ama oradaki A2 düşük (0.1–0.2); bu muhtemelen hizalanmış bir sarmal kol.
- Uzunluk medyanı 5.0 kpc (son çeyrekte 6.0) olarak kalıyor.

### 2.7 Merkezi kütle yoğunlaşması

Grafik: `runs/merkezi_kutle.png`. Ölçülen kütle, R < 1 kpc içindeki **disk parçacıklarının** kütlesi.
Sabit şişkinlik ve kara delik değişmediği için dahil değil.

| Koşu | t = 0 | 500 Myr | 1300 Myr | 1800 Myr | 2000 Myr | 1300→1800 |
|---|---|---|---|---|---|---|
| ana | 2.82 | 4.02 | 5.42 | 5.64 | 5.80 ×10⁹ M☉ | **+%4.0** |
| hafif_hale | 2.82 | 5.16 | 6.58 | 6.91 | 7.05 | +%5.0 |
| soguk | 2.82 | 6.51 | 7.10 | 7.16 | 7.11 | +%0.8 |

- Çubuklar beklendiği gibi merkeze kütle taşıyor: merkezi disk kütlesi 2 Gyr'de 2–2.5 katına çıkıyor.
- **Ama ana koşunun zayıflama döneminde artış yalnızca %4**; büyüme hızı bu dönemde ortalamanın
  yarısına düşüyor (4.6e7'ye karşı 1.1e8 M☉/100 Myr).
- **Kontrol:** hafif_hale'de merkezde daha fazla kütle birikiyor (7.0e9), ama çubuk güçlü kalıyor.
- Ana koşuda M_merkez ile A2 arasındaki −0.84'lük korelasyon, aynı dönemde biri artan biri azalan iki
  eğilimden kaynaklanıyor; nedensellik göstermiyor.
- Ayrıca burada ölçülen yoğunlaşma 1 kpc ölçeğinde. Literatürde çubuğu zayıflatan merkezi kütleler
  çok daha kompakt (≲100 pc); bu modelde böyle bir yapı yok.

**Nihai sonuç: ana koşudaki çubuğun 1300–1800 Myr arasındaki zayıflaması AÇIKLANAMADI.**

| Hipotez | Bu dönemde gözlenen | Sonuç |
|---|---|---|
| Bükülme | İç kalınlık düzgün artıyor, sıçrama yok | açıklamıyor |
| Disk ısınması | Q yalnızca %3 artıyor; asıl artış çubuk kararlıyken oluyor | açıklamıyor |
| Desen hızı kaybı | Ω_p düşmüyor (+0.56 /100 Myr) | açıklamıyor |
| Merkezi kütle | M(<1 kpc) yalnızca %4 artıyor; kontrolde daha fazla kütle var, çubuk kalıcı | açıklamıyor |

Test edilmemiş olasılıklar: sabit hale (§3), sayısal ısınma (N'ye bağlılık), çubuk–sarmal etkileşimi.
Bunlar yeni koşu gerektiriyor.

## 3. Modelin bilinen sınırları

- **Hale sabit potansiyel olarak modellendi.** Bu yüzden çubuk, dinamik sürtünmeyle haleye açısal
  momentum aktaramıyor. Canlı haleli çalışmalarda çubuklar daha güçlü ve uzun ömürlü oluyor. Bu,
  modelin bilinen bir sınırı. Ana koşudaki çubuğun zayıflaması da bu sınırdan etkilenmiş olabilir;
  bunu test etmedim (canlı hale yeni koşu gerektirir).
- **Disk ısınması kısmen sayısal.** Kalınlık N=40000'de bile N'ye duyarlıydı. Kalınlıkla ilgili sayıları
  fiziksel sonuç olarak yorumlamadan önce N=80000 ile bir kontrol önerilir.
- **hafif_hale'de disk kütlesi yeniden çözülmedi** (önceki kural gereği). Bu yüzden o modelde v_c(8.2)
  artık 229 km/s değil, daha düşük.

## 4. Açık maddeler

1. **Çubuk tanımındaki iki ek koşulun onayı** (R < 5 kpc'de arama, faz-sabit uzunluk ≥ 2 kpc).
   Değerler `analysis.py`'nin başında: `BAR_THRESHOLD`, `BAR_MIN_LENGTH`, `BAR_PHASE_TOL`, `BAR_PERSIST`.
2. **Analiz sabitleri config.py'ye taşınmalı** (CLAUDE.md §8 "sihirli sayı yok"). config'e dokunmamam
   istendiği için taşımadım.
3. **CLAUDE.md güncel değil.** Aşama 2'deki kararlar (kuadrupol, N=40000, IC'de parçacık v_c düzeltmesi)
   ve bu çubuk tanımı sözleşmeye yazılmadı.
4. **Eski biçimli koşular** (`conv_*`, `quad_*`, `diag*`, `short_200`) arayüzde görünmüyor; silinmediler.
5. **Analitik model sekmesindeki Q(R) paneli** MN diskinin M_d'sini (6.8e10) üstel Σ ile kullanıyor;
   N-cisim diski 4.82e10. Aşama 1'den kalan bir tutarsızlık.

## 5. Dosyalar

- **Kod:** `src/analysis.py`, `app.py`, `tests/test_analysis.py`
- **Grafikler (`runs/`):**
  - `a2_karsilastirma.png` — çubuk gücü ve uzunluğu
  - `cubuk_faz_profilleri.png` — A2/A0(R) ve faz(R)
  - `cubuk_zayiflama_testi.png` — bükülme ve ısınma testleri
  - `desen_hizi.png` — Ω_p(t) ve R_CR / L_çubuk
  - `merkezi_kutle.png` — M(R < 1 kpc)(t) ile çubuk gücü
  - `a2_toplu_R5.png` — eski tanım
  - `dogrulama_karsilastirma.png`
  - `<ad>/karsilastirma.png` — 0/500/1000/2000 Myr; doğrulama koşularında 0/200/400/600 Myr
  - `soguk/Q_profili.png`, `soguk/parcalanma.png`
- **Animasyonlar:** `runs/ana/animasyon.mp4`, `runs/hafif_hale/animasyon.mp4`. Sağ panel yeni çubuk
  gücünü gösteriyor.
- **Sayılar:** `runs/ozet_metrikler.json` (her koşu için `bar` bloğu dahil)
- **Yeniden üretmek için:** `.venv/bin/python -m src.analysis` (~40 s); `--no-animation` ile daha hızlı.
- **Arayüz:** `.venv/bin/streamlit run app.py --server.address 127.0.0.1`
