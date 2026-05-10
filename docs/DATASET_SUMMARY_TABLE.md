# Dataset Kullanım Özet Tablosu

## TCGA GBM 2025 Dataset Kullanım Durumu

| Dataset | Durum | Kullanıldı mı? | Örnek Sayısı | Feature Sayısı | Kullanım Amacı |
|---------|-------|---------------|--------------|----------------|-----------------|
| **metabolome_quantification.txt** | ✅ KULLANILDI | Evet | 74 | 119-136 | Metabolik reprogramming, alt tip sınıflandırma, korelasyon analizleri |
| **lipidome_positive_quantification.txt** | ✅ KULLANILDI | Evet | 74 | 336 | Lipid metabolizması, lipid-metabolite korelasyonları |
| **lipidome_negative_quantification.txt** | ✅ KULLANILDI | Evet | 74 | 250 | Lipid metabolizması (pozitif mod ile birleştirildi) |
| **clinical_patient.txt** | ✅ KULLANILDI | Evet | 73 | 23 klinik değişken | Sağkalım analizleri, nüks analizleri, prognostik modelleme |
| **mrna_seq_fpkm.txt** | ❌ KULLANILMADI | Hayır | ~100 | ~20,000 gen | README'de bahsedildi ama analize dahil edilmedi |
| **protein_quantification.txt** | ❌ KULLANILMADI | Hayır | ~100 | ~10,000 protein | README'de bahsedildi ama analize dahil edilmedi |
| **acetylprotein_quantification.txt** | ❌ KULLANILMADI | Hayır | ~100 | ~12,000 asetillenmiş protein | Mevcut ancak kullanılmadı |
| **mutations.txt** | ❌ KULLANILMADI | Hayır | 104 | 6,669 mutasyon | Mevcut ancak kullanılmadı |
| **methylation_epic.txt** | ❌ KULLANILMADI | Hayır | ~100 | ~850,000 CpG | Mevcut ancak kullanılmadı |
| **mirna_zscores.txt** | ❌ KULLANILMADI | Hayır | ~100 | 2,654 miRNA | Mevcut ancak kullanılmadı |
| **circular_rna.txt** | ❌ KULLANILMADI | Hayır | ~100 | 3,672 circRNA | Mevcut ancak kullanılmadı |

---

## Analiz İstatistikleri

### Kullanılan Veriler:
- **Toplam Örnek**: 74 GBM hastası
- **Toplam Metabolite**: 119-136
- **Toplam Lipid**: 581-582 (pozitif + negatif mod)
- **Toplam Feature**: ~700 (metabolite + lipid)
- **Sağkalım Verisi**: 73 hasta

### Analiz Sonuçları:
- **Anlamlı Korelasyonlar**: 1,676 (|r|≥0.3, FDR<0.05)
- **Alt Tipler**: 2 (KMeans, Hierarchical, NMF)
- **Diferansiyel Özellikler**: 716 (FDR<0.05, |Cohen's d|>0.5)
- **Lipid Sınıfları**: 10
- **Metabolite Kategorileri**: 6
- **ML Model Doğruluğu**: %91.3-95.7 (Random Forest)
- **Sağkalım Tahmin RMSE**: 328.82 gün

---

## GBM Kanser Bağlamında Yorum

### Kullanılan Datasetlerin Klinik Önemi:

1. **Metabolome**: 
   - GBM'nin metabolik reprogramming'ini gösterir
   - Warburg effect ve glutamin anaerobik glikoliz
   - Terapötik hedefler: Metabolik enzimler

2. **Lipidome**: 
   - Membran sentezi ve hücre bölünmesi için kritik
   - Sinyal iletimi (ceramidler, fosfolipidler)
   - Terapötik hedefler: Lipid sentez enzimleri (FASN, SCD1)

3. **Clinical Data**: 
   - Sağkalım ve prognoz
   - Nüks patternleri
   - Risk stratifikasyonu

### Kullanılmayan Datasetlerin Potansiyel Değeri:

1. **mRNA**: Metabolik yol genlerinin ekspresyonu
2. **Protein**: Metabolik enzimlerin aktivitesi
3. **Mutations**: Genotip-fenotip ilişkileri
4. **Methylation**: Epigenetik-metabolik düzenleme
5. **miRNA/circRNA**: Post-transkripsiyonel düzenleme

---

## Öneriler

### Mevcut Analiz İçin:
- ✅ Lipidomik-metabolomik entegrasyonu tamamlandı
- ✅ Alt tip sınıflandırması yapıldı
- ✅ Sağkalım modelleri geliştirildi

### Gelecek Analizler İçin:
1. **Multi-Omics Entegrasyonu**: mRNA + protein + lipid + metabolite
2. **Genotip-Fenotip**: Mutasyon + metabolik profiller
3. **Epigenetik-Metabolik**: Metilasyon + metabolik yollar
4. **Post-transkripsiyonel**: miRNA/circRNA + metabolik düzenleme

---

**Not**: Bu tablo, mevcut analizin hangi datasetleri kullandığını ve hangilerini kullanmadığını net bir şekilde göstermektedir. GBM kanser bağlamında, kullanılan datasetler metabolik reprogramming ve lipid biyolojisi açısından değerli içgörüler sağlamıştır.













