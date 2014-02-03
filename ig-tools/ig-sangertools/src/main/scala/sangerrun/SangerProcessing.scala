package sangerrun

import java.io.{FileOutputStream, IOException, FilenameFilter, File}
import common.{FileUtils, Algo}
import org.biojava.bio.program.abi.ABITrace
import alicont.AlicontFactory
import alicont.common.Scoring
import alicont.algorithms.AlgorithmType
import scala.collection.mutable.ListBuffer

class SangerProcessing(primers : Seq[(String, Boolean)]) {
  private val VL_LEADER_DEFAULT = "ATGAAATACCTGCTGCCGACCGCTGCTGCTGGTCTGCTGCTCCTCGCTGCCCAGCCGGCGATGGCTAGC"
  private val VH_LEADER_DEFAULT = "ATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCC"
  private val CONTIGS_FILENAME = "contigs.fasta"
  private val CHAINS_FILENAME_TEMPLATE = "chains-%s.fasta"
  private val _primers =
    (for (((name, reverse), idx) <- primers.zipWithIndex)
    yield name -> new {
      val index: Int = idx
      val reversed: Boolean = reverse
    }).toMap

  private case class Read(sequence: String, name: String, primer: String) {}

  def process(indir : File, outdir : File, local : Boolean = false) : Unit = {

    def touchDir(dir: File): Unit =
      if (dir.exists()) {
        if(!dir.isDirectory) throw new IOException(s"$dir is a file, not directory")
      } else if (!dir.mkdir()) throw new IOException(s"Could not create directory $dir")


    if (!indir.isDirectory) throw new IOException(s"Directory not exists: $indir")
    touchDir(outdir)

    val contigs = assemble(readSequences(indir))

    writeFasta(new File(outdir, CONTIGS_FILENAME), contigs)

    val chains = extract(contigs, local)

    var vls = new ListBuffer[(String, String)]
    var vhs = new ListBuffer[(String, String)]

    for ((vl, vh, name) <- chains) {
      vl match {
        case Some(seq) => vls += Tuple2(seq, name)
        case None => printf("Oops! VL missing o_O (%s)\n", name)
      }

      vh match {
        case Some(seq) => vhs += Tuple2(seq, name)
        case None => printf("Oops! VH missing o_O (%s)\n", name)
      }
    }

    writeFasta(new File(outdir, CHAINS_FILENAME_TEMPLATE.format("vl")), vls)
    writeFasta(new File(outdir, CHAINS_FILENAME_TEMPLATE.format("vh")), vhs)
  }


  private def writeFasta(file: File, contigs: Seq[(String, String)]) {
    val stream = new FileOutputStream(file)
    FileUtils.writeFasta(stream, contigs)
    stream.close()
  }

  private def assemble(reads: Seq[Read]) : Seq[(String, String)] = {
    val groups = reads.groupBy{ case Read(seq, name, primer) => name }
    groups.toList.map{ case (name, groupReads) =>
      val sortedReads = groupReads.sortWith((a, b) => _primers(a.primer).index < _primers(b.primer).index);
      (Assembler.getContig(sortedReads.map(r => r.sequence): _*), name)
    }
  }

  private def extract(contigs: Seq[(String, String)], local : Boolean) : Seq[(Option[String], Option[String], String)] = {
    val vlhs = contigs.map{case (seq, name) => (getVLH(seq, local), name)}
    for (((vl, vh), name) <- vlhs) yield (vl, vh, name)
  }

  private def readSequences(dir: File): Seq[Read] = {
    for (
      file <- getFiles(dir);
      split = splitFilename(file.getName)
      if split.isDefined;
      (name, primer) = split.get;
      seq = readSequence(file, rc=_primers(primer).reversed)
      if seq.isDefined
    )
      yield new Read(seq.get, name, primer)
  }

  private def splitFilename(name: String): Option[(String, String)] = { // returns (prefix, primer)
    val primerMatches = _primers.keys.map { case p => (name.indexOf(p), p) }.filter{case (i, _) => i > -1}
    if (primerMatches.isEmpty)
    None
    else {
      val (pos, primer) = primerMatches.minBy(_._1)
      Some((name.substring(0, pos), primer))
    }
  }

  private def getFiles(dir: File) : Seq[File] = {
    dir.listFiles(new FilenameFilter {
      def accept(dir: File, name: String): Boolean = name.endsWith(".ab1")
    }).sorted
  }

  private def readSequence(file : File, rc : Boolean) : Option[String] = {
    try {
      val abi = new ABITrace(file).getSequence.seqString().toUpperCase

      val r = if (!rc) abi else Algo.reverseComp(abi)
      if (!r.contains('N')) Some(r) else None
    }
    catch {
      case e : ArrayIndexOutOfBoundsException => None
    }
  }
  
  private def getVLH(seq : String, local : Boolean) : (Option[String], Option[String]) = {
    var vl_pattern : String = null
    var vh_pattern : String = null
    if (local) {
      val alicont = AlicontFactory.createSimpleAlicont(100, seq, -5,
        Scoring.loadMatrix("../../data/NUC4.4.txt"),
        AlgorithmType.SEMIGLOBAL)
      alicont.push(VL_LEADER_DEFAULT)
      val (_, vl_alignment) = alicont.alignment()
      alicont.pop()
      alicont.push(VH_LEADER_DEFAULT)
      val (_, vh_alignment) = alicont.alignment()

      vl_pattern = Algo.getAlignedPattern(vl_alignment._1, vl_alignment._2)
      vh_pattern = Algo.getAlignedPattern(vh_alignment._1, vh_alignment._2)
    }
    else {
      vl_pattern = VL_LEADER_DEFAULT
      vh_pattern = VH_LEADER_DEFAULT
    }

    val vl_lead_start = seq.indexOf(vl_pattern)
    val vh_lead_start = seq.indexOf(vh_pattern)

    var vl: Option[String] = None
    var vh: Option[String] = None

    if (vl_lead_start > -1) {
      val vl_start = vl_lead_start + vl_pattern.size
      val vl_end = if (vh_lead_start > -1) vh_lead_start else seq.length
      if (vl_start < vl_end) vl = Some(seq.substring(vl_start, vl_end))
    }

    if (vh_lead_start > -1) {
      val vh_start = vh_lead_start + vh_pattern.size
      if (vh_start < seq.length) vh = Some(seq.substring(vh_start, seq.length))
    }

    (vl, vh)
  }
}
