package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;
import android.os.Parcel;
import android.os.Parcelable;

import java.util.Date;
import java.util.List;
import java.util.Map;

/**
 * Created by itowey on 12/08/16.
 */
public class EntityDetailHistorical extends BaseObservable implements Parcelable {

    private String id;
    private String sex;
    private String age;
    private Map sire;
    private Map dam;
    private Map historicalSummary;
    private String name;
    private String owner;
    private List historicalForm;
    private Date collectionDate;

    public EntityDetailHistorical() {
    }

    protected EntityDetailHistorical(Parcel in) {
        id = in.readString();
        sex = in.readString();
        age = in.readString();
        name = in.readString();
        owner = in.readString();
        collectionDate = new Date(in.readLong());
        in.readMap(sire, this.getClass().getClassLoader());
        in.readMap(dam, this.getClass().getClassLoader());
        in.readMap(historicalSummary, this.getClass().getClassLoader());
        historicalForm = in.readArrayList(this.getClass().getClassLoader());
    }

    @Override
    public void writeToParcel(Parcel dest, int flags) {
        dest.writeString(id);
        dest.writeString(sex);
        dest.writeString(age);
        dest.writeString(name);
        dest.writeString(owner);
        dest.writeLong(collectionDate.getTime());
        dest.writeMap(sire);
        dest.writeMap(dam);
        dest.writeMap(historicalSummary);
        dest.writeList(historicalForm);
    }

    @Override
    public int describeContents() {
        return 0;
    }

    public static final Creator<EntityDetailHistorical> CREATOR = new Creator<EntityDetailHistorical>() {
        @Override
        public EntityDetailHistorical createFromParcel(Parcel in) {
            return new EntityDetailHistorical(in);
        }

        @Override
        public EntityDetailHistorical[] newArray(int size) {
            return new EntityDetailHistorical[size];
        }
    };

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSex() {
        return sex;
    }

    public void setSex(String sex) {
        this.sex = sex;
    }

    public String getAge() {
        return age;
    }

    public void setAge(String age) {
        this.age = age;
    }

    public Map getSire() {
        return sire;
    }

    public void setSire(Map sire) {
        this.sire = sire;
    }

    public Map getDam() {
        return dam;
    }

    public void setDam(Map dam) {
        this.dam = dam;
    }

    public Map getHistoricalSummary() {
        return historicalSummary;
    }

    public void setHistoricalSummary(Map historicalSummary) {
        this.historicalSummary = historicalSummary;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getOwner() {
        return owner;
    }

    public void setOwner(String owner) {
        this.owner = owner;
    }

    public List getHistoricalForm() {
        return historicalForm;
    }

    public void setHistoricalForm(List historicalForm) {
        this.historicalForm = historicalForm;
    }

    public Date getCollectionDate() {
        return collectionDate;
    }

    public void setCollectionDate(Date collectionDate) {
        this.collectionDate = collectionDate;
    }
}
