import subprocess 
import argparse
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd 
import sys
from Bio.Seq import Seq 
from ORF_nt_db import PROTEIN_TO_NUCLEOTIDE

fixed_file = open("filtered_variants.txt", "w+")

def find_aa_ref(protein,residue):
    translated_proteins = {
        "ORF1a_polyprotein": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNGFAVX",
        "ORF1ab_polyprotein_ribosomal_slippage": "RVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNNX",
        "surface_glycoprotein": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYTX",
        "ORF3a_protein": "MDLFMRIFTIGTVTLKQGEIKDATPSDFVRATATIPIQASLPFGWLIVGVALLAVFQSASKIITLKKRWQLALSKGVHFVCNLLLLFVTVYSHLLLVAAGLEAPFLYLYALVYFLQSINFVRIIMRLWLCWKCRSKNPLLYDANYFLCWHTNCYDYCIPYNSVTSSIVITSGDGTTSPISEHDYQIGGYTEKWESGVKDCVVLHSYFTSDYYQLYSTQLSTDTGVEHVTFFIYNKIVDEPEEHVQIHTIDGSSGVVNPVMEPIYDEPTTTTSVPLX",
        "envelope_protein": "MYSFVSEETGTLIVNSVLLFLAFVVFLLVTLAILTALRLCAYCCNIVNVSLVKPSFYVYSRVKNLNSSRVPDLLVX",
        "membrane_glycoprotein": "MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLWLLWPVTLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMWSFNPETNILLNVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLPKEITVATSRTLSYYKLGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALLVQX",
        "ORF6_protein": "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEIDX",
        "ORF7a_protein": "MKIILFLALITLATCELYHYQECVRGTTVLLKEPCSSGTYEGNSPFHPLADNKFALTCFSTQFAFACPDGVKHVYQLRARSVSPKLFIRQEEVQELYSPIFLIVAAIVFITLCFTLKRKTEX",
        "ORF7b": "MIELSLIDFYLCFLAFLLFLVLIMLIIFWFSLELQDHNETCHAX",
        "ORF8_protein": "MKFLVFLGIITTVAAFHQECSLQSCTQHQPYVVDDPCPIHFYSKWYIRVGARKSAPLIELCVDEAGSKSPIQYIDIGNYTVSCLPFTINCQEPKLGSLVVRCSFYEDFLEYHDVRVVLDFIX",
        "nucleocapsid_phosphoprotein": "MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQVILLNKHIDAYKTFPPTEPKKDKKKKADETQALPQRQKKQQTVTLLPAADLDDFSKQLQQSMSSADSTQAX",
        "ORF10_protein": "MGYINVFAFPFTIYSLLLCRMNSRNYIAQVDVVNFNLTX"
    }

    return translated_proteins[protein][residue-1]

def translate(seq):
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 
        'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-name', help="Provide sample name.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	sample_name = args.name

	for line in open("variants.txt"):
            dp4=line.split("DP4=")[1].split(";")[0]
            ad = line.split("AD=")[1].split(";")[0]
            #allele_ref = int(ad.split(",")[0]) + int(ad.split(",")[1])
            #allele_alt = int(ad.split(",")[2]) + int(ad.split(",")[3])
            allele_ref = int(ad.split(",")[0])
            if (int(dp4.split(",")[2]) + int(dp4.split(",")[3])) == 0:
                allele_alt = 0
            else:
                allele_alt = int(ad.split(",")[1])
            
            fixed_depth = int(dp4.split(",")[0]) + int(dp4.split(",")[1]) + int(dp4.split(",")[2]) + int(dp4.split(",")[3])

            if (allele_ref + allele_alt) > 0 and allele_alt > 0:
                if("IMF" in line):
                    af = float(line.split("IMF=")[1].split(";")[0])
                else:
                    af = allele_alt / fixed_depth
                
                type = line.split("\t")[1]

                if(af >= 0.01 and type!="synonymous SNV" and len(line.split("\t")[6])<400 and "wholegene" not in line.split("\t")[2]):
                    line_parts = line.split("\t")
                    nuc_ref = (line_parts[6])
                    nuc_alt = (line_parts[7])
                    fixed_aa_change = line_parts[2].split(":p.")[1].split(",")[0]
                    fixed_protein = line_parts[2].split(":")[1] 
                    #fixed_depth = int(allele_ref + allele_alt)
                    fixed_nuc_change = line_parts[2].split(":c.")[1].split(":")[0]
                    # if(line_parts[12].rstrip()=="-"):
                    #     mat_peptide = ""
                    #     mat_peptide_nuc_change = ""
                    #     mat_peptide_aa_change = ""
                    # else:
                    #     mat_peptide = line_parts[12].split(":")[0]
                    #     mat_peptide_nuc_change = line_parts[12].split(":")[1].rstrip().split(";")[0].strip()
                    #     mat_peptide_aa_change = line_parts[12].split(":")[1].rstrip().split(";")[1].strip()
                    #visualization.write("Sample,Position,Protein,AAChange,NucleotideChange,AlleleFreq,Depth,Type,MatPeptide,MatPeptideAAChange,MatPeptideNucChange\n"

                    nuc = fixed_nuc_change
                    if '_' in nuc:
                        nuc_num = nuc.split("_")[0]
                    elif 'del' in nuc or 'dup' in nuc:
                        nuc_num = int(nuc[0:-4])
                    elif type == "frameshift insertion" or "ins" in nuc:
                        nuc_num = int(nuc.split("_")[0])
                    else:
                        nuc_num = int(nuc[1:-1])
                    
                    nuc_num=str(nuc_num)
                                        
                    # Time to break up some indels into multiple lines
                    if(type=="nonframeshift deletion"):
                        #translated_ref = translate(nuc_ref)
                        aa_start_pos = int(fixed_aa_change.split("_")[0])
                        nuc_pos = int(line_parts[4])
                        deleted_aa = ""

                        ### these variables don't get incremented, in case we need to use the start again.
                        start_nuc_pos = nuc_pos
                        start_aa_start_pos = aa_start_pos
                        start_nuc_num = nuc_num
                        ###

                        start_split_amino_ref = find_aa_ref(fixed_protein, start_aa_start_pos)
                        
                        for codon in range(3,len(nuc_ref) + 1, 3):
                            split_nuc_ref = (nuc_ref[codon-3:codon])
                            #split_amino_ref = (translated_ref[int(codon/3)-1])
                            split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
                            deleted_aa += split_amino_ref

                            nuc_change = split_nuc_ref + str(nuc_pos) + "del"
                            #                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF              AASUB           NUCCHANGE             AAFREQ               DEPTH                                                                   TYPE
                            fixed_file.write(sample_name + "," + str(fixed_protein) + "," + str(nuc_pos) + "," + str(aa_start_pos) + "," + split_amino_ref + "," + "-" + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + str(nuc_num) + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
                            aa_start_pos += 1
                            nuc_pos +=3
                            nuc_num = int(nuc_num) + 3 
                        # if this deletion occurred in the middle of a codon, it would have resulted in a potential amino acid substitution (unless the new codon is synonymous with the ref).
                        # if not, we need to add an extra mutation row marking this.

                        base_before_del = int(start_nuc_num) - 1 # still 1-indexed, subtract 1 to get base before del
                        base_before_del_idx = base_before_del - 1 # zero indexed

                        position_in_codon = base_before_del % 3
                        if position_in_codon == 2: #the base before deletion is the last in the codon
                            continue

                        del_len = len(nuc_ref)

                        # build the amino acid with the bases following the deletion
                        codon_from_del = PROTEIN_TO_NUCLEOTIDE[fixed_protein][base_before_del_idx]

                        nuc_index = base_before_del_idx + del_len + 1
                        for i in range(2):
                            codon_from_del += PROTEIN_TO_NUCLEOTIDE[fixed_protein][nuc_index]
                            nuc_index += 1

                        # get the original amino acid contributing bases to the "right" side of the deletion
                        aa_right = find_aa_ref(fixed_protein, aa_start_pos)

                        # if different from ref, report it too.
                        aa_from_del = translate(codon_from_del.replace("U", "T"))
                        if aa_from_del != start_split_amino_ref:
                            # nuc_change = deleted_aa + aa_right + str(start_aa_start_pos) + aa_from_del + "del"
                            nuc_change = nuc_ref + str(start_nuc_pos) + "del"
                            aa_change = deleted_aa + aa_right

                            fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(start_aa_start_pos) + "," + aa_change + "," + aa_from_del + "," +  nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + "AA substitution from nonframeshift deletion"+ "," + start_nuc_num + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
                    elif(type=="nonframeshift insertion"):
                        #translated_alt = translate(nuc_alt)
                        aa_start_pos = int(fixed_aa_change.split("delins")[0][1:])
                        #split_amino_ref = fixed_aa_change[0]
                        split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
                        split_amino_alt = fixed_aa_change.split("delins")[1]
                        nuc_change = line_parts[4] + "ins" + fixed_nuc_change.split("ins")[1]
                        #                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB             NUCCHANGE             AAFREQ               DEPTH                                                                   TYPE
                        fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," +  nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1]  + "," + nuc_num + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
                    
                    elif(type == "frameshift insertion"):
                        #split_amino_ref = fixed_aa_change[0]
                        aa_start_pos = int(fixed_aa_change.split("fs")[0][1:])
                        split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
                        split_amino_alt = "fs"
                        if("dup" in nuc):
                            nuc_change = line_parts[4] + "dup" + fixed_nuc_change.split("dup")[1]
                        else:
                            nuc_change = line_parts[4] + "ins" + fixed_nuc_change.split("ins")[1]
                        #                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB               NUCCHANGE             AAFREQ               DEPTH                                                                   TYPE
                        fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
                    
                    elif(type == "frameshift deletion"):
                        #split_amino_ref = fixed_aa_change[0]
                        aa_start_pos = int(fixed_aa_change.split("fs")[0][1:])
                        split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
                        split_amino_alt = "fs"
                        nuc_change = nuc_ref + line_parts[4] + "del"

                        #                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB               NUCCHANGE             AAFREQ               DEPTH                                                                   TYPE
                        fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
                    
                    else:
                        split_amino_ref = fixed_aa_change[0]
                        split_amino_alt = fixed_aa_change[-1]
                        # Some stoplosses do not have an aasub, so check for this first...
                        if fixed_aa_change[-1].isnumeric():
                            aa_start_pos = fixed_aa_change[1:]
                        else:
                            aa_start_pos = fixed_aa_change[1:-1]
                        nuc_change = nuc_ref + line_parts[4] + nuc_alt

                        if(split_amino_ref=="X"):
                            type = "stoploss"
                            split_amino_alt = "fs"
                        else:
                            type = line_parts[1]

                        #print(sample_name)
                        #print(fixed_protein)
                        #print(line_parts[4])
                        #                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB           NUCCHANGE             AAFREQ               DEPTH                                                             TYPE
                        fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," +  nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + type + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
