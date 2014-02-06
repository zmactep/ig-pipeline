package sangerrun

import java.io.{IOException, File}
import scopt.OptionParser
import scala.util.Try

/** Created by smirnovvs on 31.01.14. */
object Main {
  private case class Config(indir: File = null, outdir: File = null, local: Boolean = false, primers: Seq[String] = Seq())

  private def getParser : scopt.OptionParser[Config] = new OptionParser[Config]("ig-sangertools") {
    head("ig-sangertools", "1.0-SNAPSHOT")
    opt[File]("indir") valueName "<dir>" action { (x, c) => c.copy(indir = x) } validate { x => if (x.isDirectory) success else failure("Not a directory")} text "input directory"
    opt[File]("outdir") valueName "<dir>" action { (x, c) => c.copy(outdir = x) } text "output directory"
    opt[Unit]("local") action { (_, c) => c.copy(local = true) } text "use local alignment to find primers in reads"
    arg[String]("<primer{+,-}>...") minOccurs 1 maxOccurs 1024 action { (x, c) => c.copy(primers = c.primers :+ x)} text "Primers in order of occurrence in plasmid. Sign indicates read direction."
    help("help") text "this message"
  }

  def main(args: Array[String]) : Unit = {
    val parser = getParser

    val config = parser.parse(args, Config()) getOrElse {
      parser.showUsage
      sys.exit()
    }

    val primers = config.primers.map(s => {
        val direction = s(s.length - 1) match {
          case '+' => false
          case '-' => true
          case _ => println("Invalid primer description format"); sys.exit()
        }
        (s.substring(0, s.length - 1), direction)
      })

    Try(new SangerProcessing(primers, config.local).process(config.indir, config.outdir)) recover {
      case e: IOException => println(e.getMessage)
      case e: Exception => println(e.printStackTrace())
    }
  }
}
