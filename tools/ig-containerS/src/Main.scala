import igcont.anno.Anno
import igcont.Container
import igcont.kmer.Counter
import igcont.trie.Trie
import scala.util.Random

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 13:10
 */

object Main{

  def trie_test() = {
    val t = new Trie()
    val rand = new Random()
    val rc = Runtime.getRuntime

    val start = System.currentTimeMillis()
    (1 until 100000).foreach(i => {
      var key = 0
      (1 to 100).foreach(j => {
        key = t.insert(key, (65+rand.nextInt(4)).toChar)
        t.setDataOf(key, rand.nextInt(100))
      })
    })

    println((System.currentTimeMillis() - start) / 1000.0)

    val m = (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024
    println(m)

    println(t.size())
  }

  def trie_iter_test() = {
    val t = new Trie()
    val sa = Array("ACGG", "GCGT", "ACGT", "GT")

    sa.foreach(s => {
      var k = 0
      s.foreach(c =>
        k = t.insert(k, c)
      )
    })

    var i = t.nextOf(0)
    while (i != 0) {
      println(i)
      i = t.nextOf(i)
    }
  }

  def kmers_test() = {
    val k = new Counter("01", '-', 3)
    val r = Array(0,1,2, 5, 3, 7, 6, 4)

    k.add(1, "0001011100", r)

    println(k.check("101", use_special = true))

    println(k.get("00-"))
    println(k.get("1-1"))
  }

  def kmers_nucleo_test(n : Int) = {
    val k = new Counter("ACGT", 'N', 5)
    val len = n - 4
    val r = 1 to len
    val rand = new Random()
    var s = ""
    (1 to n).foreach(i => s += "ACGT"(rand.nextInt(4)))
    println(s)

    k.add(1, s, r)

    println(k.get("GGNGT"))
    println(k.get("GGGGT"))
  }

  def annotation_test() = {
    val anno = new Anno(Array("Regions", "Genes", "Sites"))
    println(anno.keys())

    val rec = anno.createRecord("First", 20)
    println(rec.name())

    val i = rec.setAnnotation(0, "Regions", "FR1")
    rec.setAnnotation(1, i._1, i._2)
    rec.setAnnotation(2, "Regions", "CDR1")
    rec.setAnnotation(3, "Genes", "V")

    println(rec.annotationOf(0))
    println(rec.annotationOf(1))
    println(rec.annotationOf(2))
    println(rec.annotationOf(3))
  }

  def container_test() = {
    val cont = new Container("ACGT", 'N')

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")

    println(cont.seq(1))
    println(cont.seq("Seq1"))
  }

  def container_annotation_test() = {
    val cont = new Container("ACGT", 'N', Array("A1", "A2", "A3"), 7)

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")

    val rec = cont.record("SeqA")
    println(rec.setAnnotation(2, "A1", "42"))
    println(rec.setAnnotation(3, "A1", "42"))
    println(rec.setAnnotation(5, "A1", "38"))
    println(rec.setAnnotation(3, "A2", ":)"))

    cont.data("SeqA").foreach(tpl => {
      val (c, a) = tpl
      print(c + "  ")
      a.foreach(s => if(s._2.size > 0) print(s + " "))
      println("")
    })
  }

  def search_test() = {
    val cont = new Container("ACGT", 'N', Array(), 3)

    cont.push("GCGTGT", "Seq1")
    cont.push("AAAAAA", "SeqA")
    cont.push("AGCTTGT", "Seq2")

    val result = cont.find("CNTGT")
    println(result.size)
    println(result)
  }

  def load_test(i : Int) = {
    val cont = new Container("ACGT", 'N')
    val rand = new Random()

    val rc = Runtime.getRuntime
    val start = System.currentTimeMillis()

    for(i <- 0 until i*1000) {
      cont.push((0 until 400).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c), i.toString)
    }
    printf("%d) %.2f (%dMB)\n", i, (System.currentTimeMillis() - start) / 1000.0,
      (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024)
  }

  def warmup() = {
    def warmup_fun() = {
      val cont = new Container("ACGT", 'N')
      val rand = new Random()
      for(i <- 0 until 1000) {
        cont.push((0 until 500).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c), i.toString)
      }
    }

    (0 to 100).foreach(_ => warmup_fun())
  }

  def main(args : Array[String]) = {
//    trie_test()
//    trie_iter_test()
//    kmers_test()
//    kmers_nucleo_test(500)
//    annotation_test()
//    container_test()
//    container_annotation_test()
//    search_test()
    warmup()
    (1 to 7).foreach(i => load_test(i))
  }
}