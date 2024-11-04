12.07.2024

github'da anlatıldığı vechile indirdim ve demo'nun çalıştırılmasında bir sorun yaşamadım.
Şu an demonun ne yaptığına dair herhangi bir fikrim yok. Yalnızca üç ayrı terminalde proxy'cilik oynayan yapıların birbirini görüp konuştuklarını billiyorum.

Party.h içindeki "#include <cryptopp/osrng.h>" satırı hata veriyor gibiydi ve önerilen şekilde include path'e aşağıdaki satırın yazılması çözdü:
${workspaceFolder}/build/dependencies/cryptopp-cmake

demonun devamında söylenen FRACTIONAL_BITS'i 10'a tebdil etmeden sonra nedense make sorun çıkardı. build'i silip yeniden cmake'leyince sorun kalmadı gibi görünüyor.
Az vakit almıyor ama...

İnanılmaz saçmalıkta bir sorun var:

    uint64_t* result_of_mult = Multiply(proxy, arr1_share, arr2_share, 4); //4 tane eleman olduğunu falan da söylemek lazım.
    double* result_reconstructed = ConvertToDouble(Reconstruct(proxy, result_of_mult,4) ,4);

Bu iki satır kodda eğer "result_of_mult"u içinde "_" geçmeyen bir şekilde tesmiye edersem (misal "resul"), compiler yanlış ConvertToDouble fonksiyonunu çağırmaya çalışıyor ve kod çalışmıyor... İnanılmaz!

24.07.2024
wow, baya zaman geçmiş...

Yazdığım sigmoid fonksiyonu şu an tekil eleman alıp onun üzerinde işlem yapıyor fakat bu ideal değil. Onu vektörize etmek, yani ki bir int array almasını sağlayıp onun üzerinde koşturmak gerekiyor.

if(helper) koşulu için benim elan "sendByte()" dediğim zıkkımlar yazılmalı.
helper.cpp'ye yeni "case" eklendi. basitçe sigmoid fonksiyonunu çağırıyor.

vektörize fonksiyonlarda "param" bokunu anlamıyorum hâlâ

CrateShare metodu eğer boyut tedarik edilirse bütün bir vektörün paylarını taksim edebiliyor. 


20.09.2024
İpek'in son kodu düzgün sonuç veriyor, benimki vermiyor. Fark şuralarda olabilir:

Gradient_Descent ---- Bundan değilmiş. Acc'ın hesaplanmasında sıkıntı oluyormuş.