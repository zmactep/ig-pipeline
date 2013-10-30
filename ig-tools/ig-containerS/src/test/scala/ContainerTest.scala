import alicont.{Scoring, Alicont}
import igcont.anno.Anno
import igcont.Container

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:25
 */
object ContainerTest {
  def annotation_test() = {
    println("*** ANNOTATION TEST ***")

    val anno = new Anno(Array("Regions", "Genes", "Sites"))
    println(anno.keys)

    val rec = anno.createRecord("First", 20)
    println(rec.name)

    val i = rec.setAnnotation(0, "Regions", "FR1")
    rec.setAnnotation(1, i._1, i._2)
    rec.setAnnotation(2, "Regions", "CDR1")
    rec.setAnnotation(3, "Genes", "V")

    println(rec.annotationOf(0))
    println(rec.annotationOf(1))
    println(rec.annotationOf(2))
    println(rec.annotationOf(3))
  }

  def alignment_test() = {
    println("*** ALIGNMENT TEST ***")

    val path : String = "../../data/BLOSUM62.txt"
    val a = new Alicont("MEANLY", -5, Scoring.loadMatrix(path))
    a.push("PLE")
    a.push("ASANT")
    a.push("LY")

    println(a.target)

    val (score, alignment) = a.alignment()
    printf("%d\n%s\n%s\n", score, alignment._1, alignment._2)
  }

  def container_test() = {
    println("*** CONTAINER TEST ***")

    val cont = new Container("ACGT", 'N')

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")

    println(cont.seq(1))
    println(cont.seq("Seq1"))
  }

  def container_annotation_test() = {
    println("*** CONTAINER ANNOTATION TEST ***")

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

  def container_alignment_test() = {
    println("*** CONTAINER ALIGNMENT TEST ***")

    val path : String = "../../data/NUC4.4.txt"
    val cont = new Container("ACGT", 'N', Array("A1", "A2", "A3"), 7)

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")


    println("Top 2:")
    cont.alignment("AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -5, Scoring.loadMatrix(path), 2).foreach(align => {
      printf("Name: %s\nScore: %d (%.2f%%)\nQ: %s\nT: %s\n\n", align.name, align.score, align.similarity * 100,
        align.query, align.target)
    })

    println("More 57% similarity:")
    cont.alignment("AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -5, Scoring.loadMatrix(path), 0.57).foreach(align => {
      printf("Name: %s\nScore: %d (%.2f%%)\n\n", align.name, align.score, align.similarity * 100)
    })
  }

  def search_test() = {
    println("*** SEARCH TEST ***")

    val cont = new Container("ACGT", 'N', Array(), 3)

    cont.push("GCTGGT", "Seq1")
    cont.push("AAAAAA", "SeqA")
    cont.push("ACTTGT", "Seq2")

    val result = cont.find("NCTNGT")
    println(result.size)
    println(result)
  }
}
