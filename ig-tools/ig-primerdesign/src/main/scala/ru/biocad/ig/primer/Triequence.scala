package ru.biocad.ig.primer

import ru.biocad.ig.primer.DnaUtils.Sequence
import scala.collection.mutable.ArrayBuffer

/**
 * Created by Kos on 02.03.14.
 * ProteinTriequence = Trie + sequence (not DnaUtils.Sequence, but just string of bases).
 * Stores all possible nucleotide sequences (in quasi-trie) that represents target protein sequence.
 * Trie is "quasi", because only first two nucleotides in codon forms a trie, but alternative paths for third nucleotide
 * are grouped by nucleotide type (A, C, G or T). Ex for Leucine, we have only 4 nodes instead of 6 in the bottom:
 *     T       C
 *     |       |
 *     T    -- T
 *     |\/ /  /|
 *     |/\|  / |
 *     A  G C  T
 * ProteinTriequence supports possible synonymous substitutions in codons and possible alternative aminoacids in sequence.
 * Each nucleotide in i-th position is connected with all acceptable nucleotides for i + 1 th position
 */
sealed trait Triequence

final case class ProteinTriequence(protein: Sequence) extends Triequence {

  //index is used for O(1) access to i-th nucleotides in reverse transcribed protein sequence.
  //Otherwise we need to traverse graph from the beginning, which takes O(n)
  private val index: Vector[Vector[Node]] = {
    val aminoAcid2CodonSet: Map[String, Set[String]] = DnaUtils.aminoAcid2CodonSet
    //mutable version of index, needed only during construction
    val indexBuffer = ArrayBuffer[Vector[Node]]()

    /**
     * Builds a quasi-trie for given aminoacid and prepends to indexBuffer. Can handle alternative acids in a list,
     * but result trie depth (val layers) is always limited by one codon size (3).
     * @param acids Amino acids in j-th position of protein. Normally we have 1 aminoacid in list, but sometimes there can be an alternative for j-th position
     * @param next List of nucleotide nodes for j + 1 amino acid. Needed to connect last nucleotides of j-th acid codons with first nucleotides of j + 1 acid codons
     * @return "head" nodes of quasi-trie for j-th aminoacid. Return value used as "next" parameter for j - 1 th acid.
     */
    def processAcid(acids: List[String], next: ArrayBuffer[Node]): ArrayBuffer[Node] = {
      val layers = Vector.fill(3)(ArrayBuffer[Node]())

      def append(head: DataNode, codon: String) = {
        var current = head
        codon.zipWithIndex.foreach{case (n: Char, i: Int) =>
          //Search in node's children for 1st and 2nd nucleotide in codon, else search in "grouped" list for 3rd nucleotide.
          //See pic for Leucine above - layer for 3-rd nucleotide can have only up to 4 distinct nodes
          val nextNodeLocation = if (i < 2) current.next else layers(2)
          nextNodeLocation.find{case d: DataNode => d.nucl == n.toString case _:TerminalNode.type  => false} match {
            case Some(nextNode: DataNode) => if (i < 2) current = nextNode else current.next += nextNode
            case _ => val nextNode = DataNode(n.toString, ArrayBuffer())
              current.next += nextNode
              layers(i) += nextNode
              current = nextNode
          }
        }
      }

      val allUniqCodons: Set[String] = acids.map(a => aminoAcid2CodonSet(a)).foldLeft(Set[String]())((acc: Set[String], x: Set[String]) => acc ++ x)
      val head = DataNode("head", ArrayBuffer())
      allUniqCodons.foreach(s => append(head, s))

      //connect with previous acid in j + 1 th position
      layers(2).foreach{case n: DataNode => n.next ++= next case TerminalNode => /*impossible situation, added to suppress warning*/}
      //append 3 layers of quasi-trie to index
      layers.reverse.foreach{layer: ArrayBuffer[Node] => indexBuffer += layer.toVector}
      //return head of quasi-trie, so we can connect j - 1 th acid to current j - th acid
      layers(0)
    }

    //construct graph from last acid to first
    protein.foldRight(ArrayBuffer[Node](TerminalNode))((s: Set[String], acc: ArrayBuffer[Node]) => processAcid(s.toList, acc))
    indexBuffer.reverse.toVector
  }

  /**
   * Get nucleotide string that translates to protein
   * @param strategy given list of children selects next node to go
   * @return nucleotide string
   */
  def sample(strategy: DecisionStrategy): Option[String] = {
    def go(path: ArrayBuffer[Node], i: Int): Unit = strategy.next(path.lift(i - 1), path(i), i) map {nxt: Node => path += nxt; go(path, i + 1)}
    val dummyHead = DataNode("START", ArrayBuffer[Node]() /*to convert Vector to ArrayBuffer*/ ++ index.head)
    val path = ArrayBuffer[Node](dummyHead)
    go(path, 0)
    Some(path.tail.map{case n: DataNode => n.nucl case _: TerminalNode.type =>}.mkString)
  }

  def sampleStream(ds: DecisionStrategy): Stream[Option[String]] = Stream.cons(sample(ds), sampleStream(ds))

  /**
   * prints trie as GraphViz graph. See http://www.graphviz.org/
   * @return GraphViz-ready graph
   */
  override def toString = {
    //make GraphViz-compatible id String based on position in index
    def nodeId(node: DataNode, i: Int, j: Int): String = s"${node.nucl}_${i}_$j"

    def printClusterElement(sb: StringBuilder, i: Int, j: Int): Unit =
      if (j < index(i).size) {
        index(i)(j) match {
          case d: DataNode => {
            sb.append(s"""    node [label="${d.nucl}",shape=circle] ${nodeId(d, i, j)};\n""")
            printClusterElement(sb, i, j + 1)
          }
          case _: TerminalNode.type => ()
        }
      }

    def printLinks(sb: StringBuilder, i: Int, j: Int): Unit =
      if (j < index(i).size) {
        index(i)(j) match {
          case d: DataNode => {
            d.next.foreach{case nxt: DataNode => sb.append(s"  ${nodeId(d, i, j)} -> ${nodeId(nxt, i + 1, index(i + 1).indexWhere(that => that eq nxt/*reference equality*/))};\n") case _: TerminalNode.type =>}
            printLinks(sb, i, j + 1)
          }
          case _: TerminalNode.type => ()
        }
      }

    val sb = new StringBuilder
    sb.append(s"digraph ${protein.mkString("_").replace("(", "_").replace(")", "_").replace(", ", "")} {\n").
      append(s"""  node [fontname="verdana"];\n  fontname="Verdana";\n  rankdir="LR";\n""")

    for (j <- 0 until index.size / 3) {
      //GraphViz cluster == amino acid. Prints square around amino acid's codons
      sb.append(s"  subgraph cluster_${protein(j).mkString("_")}_$j {\n").append(s"""    label="${protein(j).mkString(", ")}";\n    graph[style=dotted];\n""")
      for (i <- 0 until 3) printClusterElement(sb, j * 3 + i, 0)
      sb.append("  }\n")
    }
    for (i <- 0 until index.size) printLinks(sb, i, 0)
    sb.append("}\n").toString()
  }
}

case object EmptyTriquence extends Triequence

object Triequence {
  def apply(seq: Option[Sequence]): Triequence = seq match {
    case Some(s) => new ProteinTriequence(s)
    case None => EmptyTriquence
  }
}

