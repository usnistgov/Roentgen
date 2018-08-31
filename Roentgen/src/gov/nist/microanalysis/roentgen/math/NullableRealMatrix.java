package gov.nist.microanalysis.roentgen.math;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.AbstractRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Creates a matrix in which one or both of the dimensions can be zero.
 * 
 * 
 * @author nicho
 *
 */
public class NullableRealMatrix extends AbstractRealMatrix {

	final int mRowDim;
	final int mColDim;
	final RealMatrix mMatrix;

	public NullableRealMatrix(int rows, int cols) {
		assert rows >= 0;
		assert cols >= 0;
		mRowDim = Math.max(0, rows);
		mColDim = Math.max(0, cols);
		if ((mRowDim > 0) && (mColDim > 0))
			mMatrix = MatrixUtils.createRealMatrix(mRowDim, mColDim);
		else
			mMatrix = null;
	}

	public NullableRealMatrix(RealMatrix rm) {
		mRowDim = rm.getRowDimension();
		mColDim = rm.getColumnDimension();
		if ((mRowDim > 0) && (mColDim > 0))
			mMatrix = rm.copy();
		else
			mMatrix = null;
	}

	public NullableRealMatrix(double[][] vals) {
		mRowDim = vals.length;
		mColDim = vals.length > 0 ? vals[0].length : 0;
		if ((mRowDim > 0) && (mColDim > 0))
			mMatrix = MatrixUtils.createRealMatrix(vals);
		else
			mMatrix = null;
	}

	private NullableRealMatrix(int rows, int cols, RealMatrix matrix) {
		mRowDim = rows;
		mColDim = cols;
		mMatrix = matrix;
	}

	@Override
	public int getRowDimension() {
		return mRowDim;
	}

	@Override
	public int getColumnDimension() {
		return mColDim;
	}

	@Override
	public RealMatrix createMatrix(int rowDimension, int columnDimension) throws NotStrictlyPositiveException {
		return new NullableRealMatrix(rowDimension, columnDimension);
	}

	@Override
	public RealMatrix copy() {
		return new NullableRealMatrix(mRowDim, mColDim, mMatrix);
	}

	@Override
	public double getEntry(int row, int column) throws OutOfRangeException {
		if (mMatrix != null)
			return mMatrix.getEntry(row, column);
		else
			throw new OutOfRangeException(0, 0, -1);
	}

	@Override
	public void setEntry(int row, int column, double value) throws OutOfRangeException {
		if (mMatrix != null)
			mMatrix.setEntry(row, column, value);
		else
			throw new OutOfRangeException(0, 0, -1);
	}
}
