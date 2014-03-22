package ru.biocad.ig.primer

import ru.biocad.ig.primer.NodeUtils._
import scala.Some

/**
 * @author kfeodorov
 **/

sealed trait DecisionStrategy {
  /**
   * Choose next nucleotide. Sequence is prepended by dummyHead node with nucl="START". It points to the list of first nucleotides in Triquence
   * @param prev previous Node. None if curr == dummyHead("START")
   * @param curr current Node. "START" for dummyHead
   * @param index current nucleotide index in protein string. 0 if if curr == dummyHead("START")
   * @return next nucleotide from "next" @param
   */
  def next(prev: Option[Node], curr: Node, index: Int): Option[Node]
}

/**
 * Simple strategy. Next nucleotide is chosen only according to it's weight
 * @param prob - map weight->nucleotide
 */
class SimpleProbabilityDecisionStrategy(prob: Map[String, Int]) extends DecisionStrategy {

  override def next(prev: Option[Node], curr: Node, index: Int): Option[Node] = {
    curr match {
      case d: DataNode =>
        val filteredProbs = prob.filter{e => d.next.map{case n: DataNode => n.nucl case _:TerminalNode.type =>}.contains(e._1)}
        if (filteredProbs.nonEmpty) findNodeByNucl(getChildren(curr), MiscUtils.sample(filteredProbs))
        else None
      case _ => None
    }
  }
}

object FreqTable {
  val total = 1581056
  val weights = Map("UUU" -> 36738, "UCU" -> 13723, "UAU" -> 26077, "UGU" -> 8732,
                    "UUC" -> 26655, "UCC" -> 14004, "UAC" -> 19204, "UGC" -> 10920,
                    "UUA" -> 22000, "UCA" -> 12367, "UAA" -> 3217, "UGA" -> 1739,
                    "UUG" -> 22187, "UCG" -> 13729, "UAG" -> 426, "UGG" -> 23955,
                    "CUU" -> 18499, "CCU" -> 11563, "CAU" -> 21480, "CGU" -> 32051,
                    "CUC" -> 17321, "CCC" -> 9154, "CAC" -> 15490, "CGC" -> 33132,
                    "CUA" -> 6256, "CCA" -> 13508, "CAA" -> 23778, "CGA" -> 6132,
                    "CUG" -> 80431, "CCG" -> 34497, "CAG" -> 46720, "CGG" -> 9988,
                    "AUU" -> 47046, "ACU" -> 14447, "AAU" -> 29426, "AGU" -> 14953,
                    "AUC" -> 38317, "ACC" -> 35978, "AAC" -> 33762, "AGC" -> 25265,
                    "AUA" -> 8607, "ACA" -> 12981, "AAA" -> 52512, "AGA" -> 4660,
                    "AUG" -> 42762, "ACG" -> 23419, "AAG" -> 16855, "AGG" -> 3041,
                    "GUU" -> 29188, "GCU" -> 24683, "GAU" -> 50786, "GGU" -> 38610,
                    "GUC" -> 23896, "GCC" -> 39754, "GAC" -> 29331, "GGC" -> 44073,
                    "GUA" -> 17622, "GCA" -> 32557, "GAA" -> 60419, "GGA" -> 14209,
                    "GUG" -> 40384, "GCG" -> 50057, "GAG" -> 27971, "GGG" -> 17812)
}

/**
 * Strategy, based on codon frequency table. Takes into account position of Node
 */
class CodonFrequencyDecisionStrategy extends DecisionStrategy{
  val equalsWeights = Map("A" -> 1000, "C" -> 1000, "G" -> 1000, "U" -> 1000)

  override def next(prev: Option[Node], curr: Node, index: Int): Option[Node] = {
    def threeNodes2SetOfCodons(start: Node): Set[String] = start match {
      case s: DataNode =>
        (for {
          child <- s.next
          grandChild <- child.asInstanceOf[DataNode].next if child.isInstanceOf[DataNode]
        } yield child match {
            case c: DataNode => grandChild match {
              case g: DataNode => s.nucl + c.nucl + g.nucl
              case _ => ""
            }
            case _ => ""
          }).filter(_.length == 3).toSet
      case _ => Set()
    }

    curr match {
      case d: DataNode =>
        if (index < 2) {
          val m = equalsWeights.filterKeys(p => d.next.map{case n: DataNode => n.nucl case _:TerminalNode.type =>}.contains(p))
          if (m.nonEmpty) findNodeByNucl(getChildren(curr), MiscUtils.sample(m))
          else None
        } else {
          index match {
            //as index == 0 for dummyHead("START"), i % 3 == 1 is the first nucl in codon
            case i if i % 3 == 1 =>
              val codons = threeNodes2SetOfCodons(curr)
              val m = FreqTable.weights.filterKeys(k => codons.contains(k))
              if (m.nonEmpty) {
                val choice = MiscUtils.sample(m)
                findNodeByNuclWithChild(getChildren(curr), choice(1).toString, choice(2).toString)
              } else None

            //second nucl in codon
            case i if i % 3 == 2 =>
              val codons: Set[String] = d.next.map{case n: DataNode => prev.get.asInstanceOf[DataNode].nucl + d.nucl + n.nucl case _:TerminalNode.type => ""}.filter(_.nonEmpty).toSet
              val m = FreqTable.weights.filterKeys(k => codons.contains(k))
              if (m.nonEmpty) {
                val choice = MiscUtils.sample(m)
                findNodeByNucl(getChildren(curr), choice(2).toString)
              } else None

            //last nucl in codon
            case i if i % 3 == 0 =>
              val codons: Set[String] = d.next.flatMap{node => threeNodes2SetOfCodons(node)}.toSet
              val m = FreqTable.weights.filterKeys(k => codons.contains(k))
              if (m.nonEmpty) {
                val choice = MiscUtils.sample(m)
                findNodeByNucl(getChildren(curr), choice(0).toString)
              } else None
          }
        }
      case _ => None
    }
  }
}