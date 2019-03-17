package gov.nist.microanalysis.roentgen.example;

import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.WriteToLaTeX;

public class CompExamples {

	public CompExamples() //
			throws ArgumentException, ParseException, IOException {
		runAnorthoclase();
		runStainless1();
		if (true) {
			runNormalization1();
			runNormalization2();
			runNormalization3();
			runNormalization4();
		}
	}

	public void runAnorthoclase() //
			throws ArgumentException, ParseException, IOException {
		final String name = "Anorthoclase";
		final Map<Composition, Number> mcn = new HashMap<>();
		final Composition sanidine = Composition.parse("K(AlSi3O8)");
		mcn.put(sanidine, new UncertainValue(0.40, 0.01));
		final Composition albite = Composition.parse("Na(AlSi3O8)");
		mcn.put(albite, new UncertainValue(0.60, 0.01));
		final Composition mix = Composition.combine(name, mcn, true);

		final Report r = new Report("Anorthoclase");
		r.addHeader("Mixture");
		r.add(albite);
		r.add(sanidine);
		r.add(mix);
		r.inBrowser(Mode.VERBOSE);

		WriteToLaTeX wtl = new WriteToLaTeX();
		try (PrintWriter pw = new PrintWriter(System.getProperty("user.home") + "\\Desktop\\Anorthoclaise.tex")) {
			pw.write("% TERSE\n");
			wtl.write(pw, mix, WriteToLaTeX.Mode.TERSE);
			pw.write("\n\n% NORMAL\n");
			wtl.write(pw, mix, WriteToLaTeX.Mode.NORMAL);
			pw.write("\n\n% VERBOSE\n");
			wtl.write(pw, mix, WriteToLaTeX.Mode.VERBOSE);
		}
	}

	public void runNormalization1() //
			throws IOException, ArgumentException {
		final Map<Element, Number> men = new HashMap<>();
		men.put(Element.Gold, new UncertainValue(0.595, 0.012));
		men.put(Element.Silver, new UncertainValue(0.402, 0.009));
		final Composition comp = Composition.massFraction("Au(60)Ag", men);
		final Report r = new Report(comp.toString());
		r.addHeader(comp.toHTML(Mode.TERSE));
		r.addSubHeader("Terse");
		r.add(comp, Mode.TERSE);
		r.addSubHeader("Normal");
		r.add(comp, Mode.NORMAL);
		r.addSubHeader("Verbose");
		r.add(comp, Mode.VERBOSE);
		r.inBrowser(Mode.VERBOSE);
		
		WriteToLaTeX wtl = new WriteToLaTeX();
		try (PrintWriter pw = new PrintWriter(System.getProperty("user.home") + "\\Desktop\\Au(60)Ag.tex")) {
			pw.write("% TERSE\n");
			wtl.write(pw, comp, WriteToLaTeX.Mode.TERSE);
			pw.write("\n\n% NORMAL\n");
			wtl.write(pw, comp, WriteToLaTeX.Mode.NORMAL);
			pw.write("\n\n% VERBOSE\n");
			wtl.write(pw, comp, WriteToLaTeX.Mode.VERBOSE);
			pw.write("\n\n% CUSTOM\n");
			List<Object> labels = new ArrayList<>();
			labels.addAll(MaterialLabel.buildMassFractionTags(comp.getMaterial()));
			// labels.addAll(Normalize.buildNormalized(MaterialLabel.buildMassFractionTags(comp.getMaterial())));
			// labels.addAll(MaterialLabel.buildAtomFractionTags(comp.getMaterial()));
			// labels.add(MaterialLabel.buildAnalyticalTotalTag(comp.getMaterial()));
			// labels.add(MaterialLabel.buildMeanAtomicNumberTag(comp.getMaterial()));
			// labels.add(MaterialLabel.buildMeanAtomicWeighTag(comp.getMaterial()));
			wtl.write(pw, comp, new BasicNumberFormat("0.000"), labels);;
		}
	}

	public void runNormalization2() //
			throws IOException, ArgumentException {
		final Map<Element, Number> men = new HashMap<>();
		men.put(Element.Gold, new UncertainValue(0.60, 0.06));
		men.put(Element.Silver, new UncertainValue(0.30, 0.03));
		men.put(Element.Copper, new UncertainValue(0.10, 0.01));
		final Composition comp = Composition.massFraction("Au(60)Ag(30)Cu(10)", men);
		final Report r = new Report(comp.toString());
		r.addHeader(comp.toHTML(Mode.TERSE));
		r.addSubHeader("Terse");
		r.add(comp, Mode.TERSE);
		r.addSubHeader("Normal");
		r.add(comp, Mode.NORMAL);
		r.addSubHeader("Verbose");
		r.add(comp, Mode.VERBOSE);
		r.inBrowser(Mode.VERBOSE);
	}

	public void runNormalization3() //
			throws IOException, ArgumentException {
		final Map<Element, Number> men = new HashMap<>();
		men.put(Element.Gold, new UncertainValue(0.60, 0.06));
		men.put(Element.Silver, new UncertainValue(0.39, 0.03));
		men.put(Element.Copper, new UncertainValue(0.01, 0.001));
		final Composition comp = Composition.massFraction("Au(60)Ag(39)Cu(1)", men);
		final Report r = new Report(comp.toString());
		r.addHeader(comp.toHTML(Mode.TERSE));
		r.addSubHeader("Terse");
		r.add(comp, Mode.TERSE);
		r.addSubHeader("Normal");
		r.add(comp, Mode.NORMAL);
		r.addSubHeader("Verbose");
		r.add(comp, Mode.VERBOSE);
		r.inBrowser(Mode.VERBOSE);
	}

	public void runNormalization4() //
			throws IOException, ArgumentException {
		final Map<Element, Number> men = new HashMap<>();
		men.put(Element.Silicon, new UncertainValue(0.532565, 0.00532565));
		men.put(Element.Oxygen, new UncertainValue(0.467435, 5.0 * 0.00467435));
		final Composition comp = Composition.massFraction("SiO2", men);
		final Report r = new Report(comp.toString());
		r.addHeader(comp.toHTML(Mode.TERSE));
		r.addSubHeader("Terse");
		r.add(comp, Mode.TERSE);
		r.addSubHeader("Normal");
		r.add(comp, Mode.NORMAL);
		r.addSubHeader("Verbose");
		r.add(comp, Mode.VERBOSE);
		r.inBrowser(Mode.VERBOSE);
	}

	public void runStainless1() //
			throws ArgumentException, IOException {
		final Material mat = new Material("SS304", Arrays.asList(Element.Chromium, Element.Nickel, Element.Iron));
		final Map<Element, Number> men = new HashMap<>();
		men.put(Element.Chromium, new UncertainValue(0.19, 0.01));
		men.put(Element.Nickel, new UncertainValue(0.0925, 0.0125));
		final Report r = new Report("Element-by-Difference");
		WriteToLaTeX wtl = new WriteToLaTeX();
		try {
			{
				final Composition comp = Composition.elementByDifference(mat, Element.Iron, men, Collections.emptyMap());
				r.add(comp, Mode.NORMAL);
				r.add(comp, Mode.VERBOSE);
				try (PrintWriter pw = new PrintWriter(System.getProperty("user.home") + "\\Desktop\\SS304.tex")) {
					pw.write("% TERSE\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.TERSE);
					pw.write("\n\n% NORMAL\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.NORMAL);
					pw.write("\n\n% VERBOSE\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.VERBOSE);
					pw.write("\n\n% CUSTOM\n");
					List<Object> labels = new ArrayList<>();
					labels.addAll(MaterialLabel.buildMassFractionTags(comp.getMaterial()));
					// labels.addAll(Normalize.buildNormalized(MaterialLabel.buildMassFractionTags(comp.getMaterial())));
					// labels.addAll(MaterialLabel.buildAtomFractionTags(comp.getMaterial()));
					// labels.add(MaterialLabel.buildAnalyticalTotalTag(comp.getMaterial()));
					// labels.add(MaterialLabel.buildMeanAtomicNumberTag(comp.getMaterial()));
					// labels.add(MaterialLabel.buildMeanAtomicWeighTag(comp.getMaterial()));
					wtl.write(pw, comp, new BasicNumberFormat("0.000"), labels);;
				}
			}
			{
				final Composition comp = Composition.elementByDifference(mat, Element.Iron, men, Collections.emptyMap());
				comp.setCalculator(new UncertainValuesCalculator.FiniteDifference(comp.getInputValues().mapMultiply(0.001)));
				r.add(comp, Mode.NORMAL);
				r.add(comp, Mode.VERBOSE);
			}
			{
				r.addHeader("Naive - mass fraction model");
				men.put(Element.Iron, new UncertainValue(0.7175,0.0160));
				final Composition comp = Composition.massFraction(mat, men, Collections.emptyList());
				r.add(comp, Mode.NORMAL);
				r.add(comp, Mode.VERBOSE);
				try (PrintWriter pw = new PrintWriter(System.getProperty("user.home") + "\\Desktop\\SS304_naive.tex")) {
					pw.write("% TERSE\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.TERSE);
					pw.write("\n\n% NORMAL\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.NORMAL);
					pw.write("\n\n% VERBOSE\n");
					wtl.write(pw, comp, WriteToLaTeX.Mode.VERBOSE);
					pw.write("\n\n% CUSTOM\n");
					List<Object> labels = new ArrayList<>();
					labels.addAll(MaterialLabel.buildMassFractionTags(comp.getMaterial()));
					// labels.addAll(Normalize.buildNormalized(MaterialLabel.buildMassFractionTags(comp.getMaterial())));
					// labels.addAll(MaterialLabel.buildAtomFractionTags(comp.getMaterial()));
					// labels.add(MaterialLabel.buildAnalyticalTotalTag(comp.getMaterial()));
					// labels.add(MaterialLabel.buildMeanAtomicNumberTag(comp.getMaterial()));
					// labels.add(MaterialLabel.buildMeanAtomicWeighTag(comp.getMaterial()));
					wtl.write(pw, comp, new BasicNumberFormat("0.000"), labels);;
				}
			}
			
		} finally
		{
			r.inBrowser(Mode.VERBOSE);
		}
	}

	public static void main(final String[] args) {
		try {
			new CompExamples();
		} catch (ArgumentException | ParseException | IOException e) {
			e.printStackTrace();
		}
	}

}
