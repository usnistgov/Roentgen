package gov.nist.microanalysis.roentgen.utility;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.models.Normalize;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;

/**
 * Helpful methods for outputing data in LaTeX friendly formats.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class WriteToLaTeX {

	private BasicNumberFormat mMassFractionFormat;
	private BasicNumberFormat mAtomFractionFormat;
	private BasicNumberFormat mAtomWeightFormat;
	private BasicNumberFormat mAtomNumberFormat;
	private BasicNumberFormat mGenericFormat;

	/**
	 *
	 */
	public WriteToLaTeX() {
		setMassFractionFormat(new BasicNumberFormat("0.0000"));
		setAtomFractionFormat(new BasicNumberFormat("0.0000"));
		setAtomWeightFormat(new BasicNumberFormat("0.00"));
		setAtomNumberFormat(new BasicNumberFormat("0.00"));
		setGenericFormat(new BasicNumberFormat("0.00E0"));
	}

	public enum Mode {
		TERSE, NORMAL, VERBOSE
	};

	public void write(final PrintWriter wr, final UncertainValues uvs, final BasicNumberFormat bnf,
			final List<Object> labels) {
		wr.write("\\begin{equation}\n");
		wr.write("\\begin{matrix}\n");
		wr.write(" \\begin{bmatrix}\n");
		wr.write(" Quantity & Value \\\\\n ");
		for (int r = 0; r < labels.size(); ++r) {
			wr.write("\\mathit{" + labels.get(r).toString() + "}");
			wr.write(" & ");
			wr.write(bnf.formatLaTeX(uvs.getEntry(labels.get(r))));
			if (r != labels.size() - 1)
				wr.write("\\\\");
			wr.write("\n ");
		}
		wr.write("\\end{bmatrix}\n");
		wr.write(" & \\pm & \n");
		wr.write(" \\begin{bmatrix}\n ");
		for (int c = 0; c < labels.size(); ++c) {
			wr.write("\\mathit{" + labels.get(c).toString() + "}");
			wr.write(c < labels.size() - 1 ? " & " : "\\\\");
		}
		wr.write("\n ");
		final BasicNumberFormat cfmt = new BasicNumberFormat("0.00");
		for (int r = 0; r < labels.size(); ++r) {
			final int rIdx = uvs.indexOf(labels.get(r));
			assert rIdx >= 0;
			for (int c = 0; c < labels.size(); ++c) {
				final int cIdx = uvs.indexOf(labels.get(c));
				assert cIdx >= 0;
				if (r == c)
					wr.write("(" + bnf.format(uvs.getUncertainty(cIdx)) + ")^2");
				else {
					final double corrCoeff = uvs.getCorrelationCoefficient(rIdx, cIdx);
					if (Math.abs(corrCoeff) < 1.0e-4) {
						wr.write("-");
					} else {
						wr.write("\\num{" + cfmt.format(corrCoeff) + "}");
						wr.write("\\cdot \\sigma_R \\sigma_C");
					}
				}
				wr.write(c < labels.size() - 1 ? " & " : "");
			}
			if (r != labels.size() - 1)
				wr.write("\\\\");
			wr.write("\n ");
		}
		wr.write(" \\end{bmatrix}\n");
		wr.write("\\end{matrix}\n");
		wr.write("\\end{equation}\n");
		wr.flush();
	}

	public void write(final PrintWriter wr, final Composition comp, final Mode mode) {
		switch (mode) {
		case TERSE: {
			wr.write(stripName(comp));
			wr.write(" = [");
			boolean first = true;
			for (final Element elm : comp.getElementSet()) {
				if (!first)
					wr.write(", ");
				final UncertainValue mf = comp.getMassFraction(elm);
				wr.write("(" + elm.getAbbrev() + ", " + mMassFractionFormat.formatLaTeX(mf) + ")");
				first = false;
			}
			wr.write("\\ by\\ mass.]");
		}
			break;
		case NORMAL: {
			wr.write("\\begin{table}[h]\n");
			wr.write("\\begin{center}\n");
			wr.write("\\begin{tabular}{lccccc}\n");
			wr.write(
					"\\textit{Element} & \\textit{Z} & \\textit{A} & \\textit{Mass Fraction} & \\textit{Normalized Mass} & \\textit{Atom Fraction} \\\\\n");
			for (final Element elm : comp.getElementSet()) {
				final UncertainValue mf = comp.getMassFraction(elm);
				final UncertainValue nf = comp.getNomalizedMassFraction(elm);
				final UncertainValue af = comp.getAtomFraction(elm);
				final UncertainValue a = comp.getAtomicWeight(elm);
				wr.write(elm.toString());
				wr.write(" & ");
				wr.write("\\num{" + Integer.toString(elm.getAtomicNumber()) + "}");
				wr.write(" & ");
				wr.write(mAtomWeightFormat.formatLaTeX(a));
				wr.write(" & ");
				wr.write(mMassFractionFormat.formatLaTeX(mf));
				wr.write(" & ");
				wr.write(mMassFractionFormat.formatLaTeX(nf));
				wr.write(" & ");
				wr.write(mAtomFractionFormat.formatLaTeX(af));
				wr.write("\\\\\n");
			}
			wr.write("\\textit{Mean/Total} & ");
			wr.write(mAtomNumberFormat.formatLaTeX(comp.getMeanAtomicNumber()));
			wr.write(" & ");
			wr.write(mAtomWeightFormat.formatLaTeX(comp.getMeanAtomicWeight()));
			wr.write(" & ");
			wr.write(mMassFractionFormat.formatLaTeX(comp.getAnalyticalTotal()));
			wr.write(" & ");
			wr.write(mMassFractionFormat.formatLaTeX(UncertainValue.ONE));
			wr.write(" & ");
			wr.write(mAtomFractionFormat.formatLaTeX(UncertainValue.ONE));
			wr.write("\n");
			wr.write("\\end{tabular}\n");
			wr.write("\\end{center}\n");
			wr.write("\\label{tbl:comp_" + stripName(comp).replace(" ", "-") + "}\n");
			wr.write("\\caption{Composition of " + stripName(comp) + "}\n");
			wr.write("\\end{table}");
		}
			break;
		case VERBOSE:
			final List<Object> labels = new ArrayList<>();
			labels.addAll(MaterialLabel.buildMassFractionTags(comp.getMaterial()));
			labels.addAll(Normalize.buildNormalized(MaterialLabel.buildMassFractionTags(comp.getMaterial())));
			labels.addAll(MaterialLabel.buildAtomFractionTags(comp.getMaterial()));
			labels.add(MaterialLabel.buildAnalyticalTotalTag(comp.getMaterial()));
			labels.add(MaterialLabel.buildMeanAtomicNumberTag(comp.getMaterial()));
			labels.add(MaterialLabel.buildMeanAtomicWeighTag(comp.getMaterial()));
			write(wr, comp, mMassFractionFormat, labels);
			break;
		}
	}

	private String stripName(final Composition comp) {
		return stripName(comp.toString());
	}

	private String stripName(final String label) {
		return label.substring(label.indexOf('[')+1, label.length() - 1);
	}

	public BasicNumberFormat getMassFractionFormat() {
		return mMassFractionFormat;
	}

	public void setMassFractionFormat(final BasicNumberFormat massFractionFormat) {
		this.mMassFractionFormat = massFractionFormat;
	}

	public BasicNumberFormat getAtomFractionFormat() {
		return mAtomFractionFormat;
	}

	public void setAtomFractionFormat(final BasicNumberFormat atomFractionFormat) {
		this.mAtomFractionFormat = atomFractionFormat;
	}

	public BasicNumberFormat getAtomWeightFormat() {
		return mAtomWeightFormat;
	}

	public void setAtomWeightFormat(final BasicNumberFormat atomWeightFormat) {
		this.mAtomWeightFormat = atomWeightFormat;
	}

	public BasicNumberFormat getAtomNumberFormat() {
		return mAtomNumberFormat;
	}

	public void setAtomNumberFormat(final BasicNumberFormat atomNumberFormat) {
		this.mAtomNumberFormat = atomNumberFormat;
	}

	public BasicNumberFormat getGenericFormat() {
		return mGenericFormat;
	}

	public void setGenericFormat(final BasicNumberFormat genericFormat) {
		this.mGenericFormat = genericFormat;
	}
}
