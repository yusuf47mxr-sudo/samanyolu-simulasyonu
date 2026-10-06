# Samanyolu Dinamiği

Samanyolu'nun dinamiğini simüle eden ve fiziksel parametreleri canlı değiştirilebilen bir arayüzle
sunan bir Python projesi. Bilimsel mesajı tek cümle: **galaksiyi bir arada tutan şey görünür madde değil.**

![3B koşu sonuçları: hafif haleli galakside çubuk, t = 2000 Myr](docs/ekran_3b.png)

Proje iki aşamada kuruldu:

1. **Test parçacık modeli:** Disk (Miyamoto–Nagai), şişkinlik (Hernquist), karanlık madde halesi (NFW)
   ve merkezi kara delikten oluşan analitik potansiyelde dönme eğrisi ve yıldız yörüngeleri.
2. **N-cisim modeli:** 40 000 parçacıklı, öz-yerçekimli bir yıldız diski. Hale sabit potansiyeldir.
   Çubuk ve sarmal kollar kendiliğinden oluşur.

## Ana bulgular

**Dönme eğrisi gözlemle uyumlu.** Varsayılan modelde Güneş'in konumunda (8.2 kpc) dönme hızı
**228.7 km/s**, yörünge periyodu **220 Myr**. Eğri 5–20 kpc arasında 212–232 km/s bandında düz kalıyor.

**Karanlık madde olmadan galaksi bir arada durmuyor.** Hale kapatılınca 20 kpc'deki dönme hızı
212 km/s'den **123 km/s**'ye düşüyor (Kepler düşüşü). Arayüzün "Canlı" sekmesinde haleyi kapatınca dış
disk birkaç yüz Myr içinde dağılıyor.

**Merkezdeki kara delik dönme eğrisini etkilemiyor.** Kütlesini 100 katına çıkarmak Güneş'in konumundaki
hızı yalnızca 0.5 km/s değiştiriyor.

**Hale hafifleyince çubuk daha erken ve daha güçlü oluşuyor.** Disk kütlesi aynı kalıp hale kütlesi
yarıya indirildiğinde:

| | Normal hale | Hafif hale (M_vir / 2) |
|---|---|---|
| Çubuk oluşumu | 390 Myr | **150 Myr** |
| En yüksek çubuk gücü (A2/A0) | 0.38 | **0.50** |
| Çubuk uzunluğu (medyan) | 3.5 kpc | **5.0 kpc** |
| 2000 Myr'de çubuk | büyük ölçüde sönmüş | **kalıcı** |

Hafif haleli çubuğun eşdönme yarıçapı / çubuk uzunluğu oranı **1.11**. Bu, gerçek gökadalarda
gözlenen 1.0–1.4 aralığının içinde.

**Kararsız disk kendini düzenliyor.** Toomre Q = 0.8 ile başlatılan disk ilk 200 Myr'de topaklara
bölünüyor, sonra kendini ısıtıyor; Q 800 Myr'de 2.3'e çıkıyor.

Ayrıntılı sonuçlar, yöntemler ve açık sorular: [runs/OZET.md](runs/OZET.md).

![Analitik model: dönme eğrisi ve bileşenleri](docs/ekran_analitik.png)

## Doğruluk

- **Kabul testleri:** 54 otomatik test var (`tests/`), aralarında projenin kabul kriterleri T1–T9 da var.
  Kriterler [CLAUDE.md](CLAUDE.md)'de.
- **Zaman entegrasyonu:** leapfrog (kick–drift–kick). 2000 Myr'lik N-cisim koşularında enerji sapması
  3e-5 – 8e-5.
- **Kuvvet hesabı:** kuadrupol momentli Barnes–Hut ağacı (numba). RMS kuvvet hatası ~%0.1. Ağaçla elde edilen
  sonuçlar doğrudan N² toplamayla karşılaştırıldı; fark rastgele başlangıç gürültüsünden küçük.

## Çalıştırma

Python 3.12 ile denendi.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

.venv/bin/python -m pytest tests                               # testler
.venv/bin/streamlit run app.py --server.address 127.0.0.1      # arayüz
```

Arayüzde üç sekme var:

- **Analitik model:** kaydıraklarla dönme eğrisi, kütleler, Toomre Q ve örnek bir yörünge.
- **Canlı:** 5000 test parçacığı gerçek zamanlı ve 3B olarak dönüyor; kaydıraklar anında etkiliyor.
- **Koşu sonuçları:** önceden koşulmuş N-cisim simülasyonlarının 3B görünümü, çubuk gücü ve zaman kaydırağı.

**"Koşu sonuçları" sekmesi için simülasyon verisi gerekir.** Tam snapshot'lar ~930 MB olduğu için depoya
konmadı. Depoda bunun yerine küçültülmüş bir kopya var (`deploy/runs/`, ~9 MB): üç koşu, her 8. snapshot,
5000 parçacık. Çubuk gücü ve enerji eğrileri yine tam çözünürlüklü koşudan geliyor. Bu veriyle çalıştırmak
için (bulutta da bu giriş noktası kullanılır; canlı modda 2000 parçacık):

```bash
.venv/bin/streamlit run deploy/app.py
```

Tam koşuları yeniden üretmek için:

```bash
.venv/bin/python -m src.nbody --name ana --t_end 2000 --N 40000          # ~25 dk
.venv/bin/python -m src.nbody --name soguk --t_end 2000 --N 40000 --Q 0.8
.venv/bin/python -m src.analysis                                         # grafikler ve animasyonlar
.venv/bin/python deploy/make_deploy_data.py                              # küçültülmüş kopya
```

Hafif haleli koşu komut satırından değil, Python'dan başlatılır:
`src.nbody.run("hafif_hale", GalaxyParams(halo=HaloParams(M_vir=5e11)), num=NumericalParams(N=40000), t_end=2000)`.

## Dosya yapısı

```
app.py            Streamlit arayüzü
src/              fizik ve sayısal motor (birimler: kpc, Myr, M☉)
  potentials.py   analitik potansiyeller
  rotation.py     dönme eğrisi, kappa, Toomre Q
  orbits.py       leapfrog
  ics.py          N-cisim başlangıç koşulları
  nbody.py        Barnes–Hut ağacı (kuadrupol), N-cisim koşusu
  analysis.py     çubuk gücü, desen hızı, animasyon
  config.py       tüm varsayılan parametreler
ui/galaxy3d/      3B görünüm bileşeni
tests/            otomatik testler
runs/             simülasyon çıktıları (yalnızca özet ve seçili grafikler depoda)
deploy/           bulut dağıtımı: küçültülmüş veri ve giriş noktası
docs/             ekran görüntüleri
```

## Bilinen sınırlar

- **Hale sabit potansiyel.** Bu yüzden çubuk, dinamik sürtünmeyle haleye açısal momentum
  aktaramıyor. Canlı haleli simülasyonlarda çubuklar daha güçlü ve uzun ömürlü oluyor.
- **Normal haleli koşuda çubuğun zayıflaması açıklanamadı.** Bükülme, disk ısınması, desen hızı kaybı
  ve merkezi kütle birikmesi test edildi; hiçbiri zamanlamayı açıklamıyor.
- **Disk kalınlığı parçacık sayısına duyarlı.** 40 000 parçacıkta bile tam yakınsamış değil.

## Lisans

[MIT](LICENSE)
